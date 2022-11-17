setwd("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Data/adding Targ to genomewide/mutants/")
library(data.table)
rm(list=ls()); gc()

files <- as.list(setdiff(list.files(pattern = "tsv"), list.dirs(recursive = FALSE, full.names = FALSE)))
#----------
compileData <- function(fileName) {
  tissue <-
    gsub(
      pattern = "_tsv",
      replacement = "",
      x = paste0(
        unlist(
          stringi::stri_extract_all_regex(str = gsub("del", "DEL", fileName , fixed = T), pattern = "[[:lower:]]+[_]*"),
          use.names = F
        ),
        collapse = ""
      )
    )
  
  mutID <-
    gsub(
      pattern = "^[_]+",
      replacement = "",
      paste0(unlist(
        stringi::stri_extract_all_regex(str = gsub("del", "DEL", fileName , fixed = T), pattern = "[^a-z.]+")
      ), collapse = ""),
      fixed = F
    )
  # reading file
  myDT <- data.table::fread(file = fileName ,
        sep = "\t",
        skip = 0,
        header = T,
        stringsAsFactors = F,
        showProgress = T,
        nThread = 1)

  # Keeping rows where there are at least 10 samples in either "threshold" samples in 
  # gw or full i.e., gw+targ data. AND where there is at least 1 sample analyzed in the gw data.
  myDT <- myDT[which((!is.na(myDT[, 2])) & 
                 (myDT[, 3] > 10 | myDT[, 6] > 10)), ]
  
  colnames(myDT)[3] <- "gw_mutant_N"
  colnames(myDT)[6] <- "full_mutant_N"
  
  myDT$mutID <- mutID
  myDT$tissue <- tissue
  return(myDT)
}

myDT <- compileData(files[[1]])
for (fileName in files[2:length(files)]){
  myDT <- rbind(myDT,compileData(fileName = fileName))
}

mutInfo <- stringi::stri_split_fixed(str = myDT$mutID,pattern = "_",simplify = T)
myDT$Gene <- mutInfo[,1]
myDT$Mutation <- gsub("_$","",paste0(mutInfo[,2],"_",mutInfo[,3]))


source("/Users/deepankar/OneDrive - O365 Turun yliopisto/Git/Gitlab.DC/Utilities/SignedFoldChange.R")
myDT[is.na(myDT$gw_percentage), gw_percentage := 0]
myDT[, FC := FoldChange(gw_percentage, full_percentage)]
myDT[, DIFF := full_percentage-gw_percentage]
myDT[is.na(myDT$gw_mutant_N), gw_mutant_N := 0]

# Reading color scheme
colors_data <- data.table::fread(file = "../color_vector.tsv")
colors_data[,tissue:=gsub(" ","_",tissue)]
colors_data <- colors_data[colors_data$tissue%in%myDT$tissue]
colors <- colors_data$color
names(colors) <- colors_data$tissue

library(ggplot2)
library(patchwork)
source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type = "L",legend = F)

#---------
main_plot <- ggplot(myDT,
       aes(
         x = FC,
         y = DIFF,
         fill = tissue,
         size = full_mutant_N,
         alpha = gw_mutant_N
       )) +
  geom_point(shape=21, stroke = 0.2, color="#000000") +
  xlab("Fold change (+/-)\n(Full vs Genome-wide)")+
  ylab("Difference in MAF (%)\n(Full - Genome-wide)")+
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = colors)+
  # scale_color_viridis_d(option = "turbo","Primary-Site") +
  scale_alpha("N (gw) / opacity",range = c(0.15, 0.8))+
  scale_size("N (full) / size", range=c(1,8),trans = "log10")+
  coord_flip(clip = "off")

scatterPlot=main_plot+customtheme
legend=ggpubr::as_ggplot(ggpubr::get_legend(main_plot))

export <-  scatterPlot + 
  legend + 
  plot_layout(nrow = 1, ncol=2,
              widths =  unit(c(8, 1), c('cm', 'null')),
              heights = unit(c(8, 1), c('cm', 'null')))

ggsave(
  plot = export,
  filename = "../mostDeviantHistologies.pdf",
  device = cairo_pdf,
  width = 24,
  height = 24,
  units = "cm"
)

#------------
# 
# main_plot <- ggplot(myDT,
#                     aes(
#                       x = gw_mutant_N,
#                       y = DIFF,
#                       fill = tissue,
#                       size = gw_N)) +
#   geom_point(shape=21, 
#              color="#000000", 
#              stroke = 0.2,
#              alpha = 1) +
#   xlab("Number samples with mutation")+
#   ylab("Difference in MAF (%)\n(Full - Genome-wide)")+
#   geom_hline(yintercept = 0) +
#   scale_fill_manual(values = colors)+
#   scale_x_continuous(breaks = seq(0,1000, by  = 100), expand = c(0,0))+
#   scale_size("Scale of study (N)", range=c(1,10))+
#   coord_flip(clip = "off")
# 
# library(patchwork)
# scatterPlot=main_plot+customtheme
# legend=ggpubr::as_ggplot(ggpubr::get_legend(main_plot))
# 
# export <-  scatterPlot + 
#   legend + 
#   plot_layout(nrow = 1, ncol=2,
#               widths =  unit(c(12, 1), c('cm', 'null')),
#               heights = unit(c(12, 1), c('cm', 'null')))
# 
# ggsave(
#   plot = export,
#   filename = "../mostDeviantHistologies2.pdf",
#   device = cairo_pdf,
#   width = 24,
#   height = 24,
#   units = "cm"
# )

#-------------
main_plot <- ggplot(myDT,
                    aes(
                      x = gw_N,
                      y = DIFF,
                      fill = tissue,
                      size = gw_mutant_N)) +
  geom_point(shape=21, 
             color="#000000", 
             stroke = 0.2,
             alpha = 1) +
  xlab("Samples analyzed (N)\n(with WXS/WGS)")+
  ylab("Difference in MAF (%)\n(Full - Genome-wide)")+
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = colors)+
  scale_size("Number of mutants", range=c(1,8))+
  scale_x_log10()+
  coord_flip(clip = "off")+
  annotation_logticks(sides = "l", outside = T)

library(patchwork)
scatterPlot=main_plot+customtheme
legend=ggpubr::as_ggplot(ggpubr::get_legend(main_plot))

export <-  scatterPlot + 
  legend + 
  plot_layout(nrow = 1, ncol=2,
              widths =  unit(c(8, 1), c('cm', 'null')),
              heights = unit(c(8, 1), c('cm', 'null')))

ggsave(
  plot = export,
  filename = "../mostDeviantHistologies3.pdf",
  device = cairo_pdf,
  width = 24,
  height = 24,
  units = "cm"
)

#-----------
# main_plot <- ggplot(myDT,
#                     aes(
#                       x = gw_mutant_N,
#                       y = FC,
#                       fill = tissue,
#                       size = gw_N)) +
#   geom_point(shape=21, 
#              color="#000000", 
#              stroke = 0.2,
#              alpha = 1) +
#   xlab("Number samples with mutation")+
#   ylab("Fold change (+/-)\n(Full vs Genome-wide)")+
#   geom_hline(yintercept = 0) +
#   scale_fill_manual(values = colors)+
#   scale_x_continuous(breaks = seq(0,1000, by  = 100), expand = c(0,0))+
#   scale_size("Scale of study (N)", range=c(1,10))+
#   coord_flip(clip = "off")
# 
# library(patchwork)
# scatterPlot=main_plot+customtheme
# legend=ggpubr::as_ggplot(ggpubr::get_legend(main_plot))
# 
# export <-  scatterPlot + 
#   legend + 
#   plot_layout(nrow = 1, ncol=2,
#               widths =  unit(c(12, 1), c('cm', 'null')),
#               heights = unit(c(12, 1), c('cm', 'null')))
# 
# ggsave(
#   plot = export,
#   filename = "../mostDeviantHistologies4.pdf",
#   device = cairo_pdf,
#   width = 24,
#   height = 24,
#   units = "cm"
# )

#-------------
main_plot <- ggplot(myDT,
                    aes(
                      x = gw_N,
                      y = FC,
                      fill = tissue,
                      size = gw_mutant_N)) +
  geom_point(shape=21, 
             color="#000000", 
             stroke = 0.2,
             alpha = 1) +
  xlab("Samples analyzed (N)\n(with WXS/WGS)")+
  ylab("Fold change (+/-)\n(Full vs Genome-wide)")+
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = colors)+
  scale_x_continuous(breaks = seq(0,3000, by  = 500), expand = c(0,0))+
  scale_size("Number of mutants", range=c(1,8))+
  coord_flip(clip = "off")

library(patchwork)
scatterPlot=main_plot+customtheme
legend=ggpubr::as_ggplot(ggpubr::get_legend(main_plot))

export <-  scatterPlot + 
  legend + 
  plot_layout(nrow = 1, ncol=2,
              widths =  unit(c(8, 1), c('cm', 'null')),
              heights = unit(c(8, 1), c('cm', 'null')))

ggsave(
  plot = export,
  filename = "../mostDeviantHistologies5.pdf",
  device = cairo_pdf,
  width = 24,
  height = 24,
  units = "cm"
)

main_plot <- ggplot(myDT,
                    aes(
                      x = gw_N,
                      y = FC,
                      fill = tissue,
                      size = gw_mutant_N)) +
  geom_point(shape=21, 
             color="#000000", 
             stroke = 0.2,
             alpha = 1) +
  xlab("Samples analyzed (N)\n(with WXS/WGS)")+
  ylab("Fold change (+/-)\n(Full vs Genome-wide)")+
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = colors)+
  # scale_x_continuous(breaks = seq(0,3000, by  = 500), expand = c(0,0))+
  scale_x_log10()+
  scale_size("Number of mutants", range=c(1,8))+
  coord_flip(clip = "off")+
  annotation_logticks(sides = "l", outside = T)

library(patchwork)
scatterPlot=main_plot+customtheme
legend=ggpubr::as_ggplot(ggpubr::get_legend(main_plot))

export <-  scatterPlot + 
  legend + 
  plot_layout(nrow = 1, ncol=2,
              widths =  unit(c(8, 1), c('cm', 'null')),
              heights = unit(c(8, 1), c('cm', 'null')))

ggsave(
  plot = export,
  filename = "../mostDeviantHistologies5_log.pdf",
  device = cairo_pdf,
  width = 24,
  height = 24,
  units = "cm"
)

# Writing data
write.table(
  x = myDT,
  "../mostDeviantHistologies.tsv",
  sep = "\t",
  row.names = F,
  col.names = ,
  quote = F
)



myDT[,Mutant := paste0(Gene," ",Mutation)]

library(plotly)
p <- plot_ly(
  data = myDT,
  size = ~gw_mutant_N,
  sizes = c(5, 40),
  x =  ~FC,
  y =  ~gw_N,
  color = ~tissue,
  colors = colors,
  mode = "text",
  text =  ~ paste("Mutant: ", paste0(Mutant," / ",gsub("_"," ",tissue)," - ", gsub("_"," ",Histology)),
                  "\n(Full): ",
                  paste0(formatC(signif(full_percentage, digits = 3), digits = 3, format = "fg"),
                         " % (n = ",
                         full_mutant_N,
                         ")"),
                  "\n (GW): ", 
                  paste0(formatC(signif(gw_percentage, digits = 3), digits = 3, format = "fg"), 
                         " % (n = ", 
                         gw_mutant_N, 
                         ")"))) %>% # write counts
  add_markers(marker=list(opacity = 0.75, sizemode = 'diameter')) %>%
  layout(
    title = paste0("<b>Most deviant histologies</b>"),
    xaxis = list(title = paste0("<b>MAF in genome-wide screens (%)</b>")),
    yaxis = list(title = paste0("<b>MAF in combined screen data (%)</b>"))
  )

# p %>% add_trace(x = ~full, y = fitted(lm(selection$gw~selection$full), mode = "lines"))

# fv <- fitted(lm(selection$gw~selection$full))

# hide_legend(p)

htmlwidgets::saveWidget(
  widget = hide_legend(p),
  file = "../mostDeviantHistologies3.html",
  selfcontained = T,
  libdir = NULL)


# plotDT <- myDT[(FC >=5 & FC < Inf)|(FC <= -5 & FC > -Inf),]

myDT[, FC := full_percentage/gw_percentage]
main_plot <- ggplot(myDT,
                    aes(
                      x = gw_N,
                      y = FC,
                      fill = tissue,
                      size = gw_mutant_N)) +
  geom_point(shape=21, 
             color="#000000", 
             stroke = 0.2,
             alpha = 1) +
  xlab("Samples analyzed (N)\n(with WXS/WGS)")+
  ylab("Fold change (+/-)\n(Full vs Genome-wide)")+
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = colors)+
  scale_x_continuous(breaks = seq(0,3000, by  = 500), expand = c(0,0))+
  geom_hline(yintercept = c(1))+
  scale_y_log10(breaks=c(1/16,1/8,1/4,1/2,1,2,4,8,16), labels=c("1/16","1/8","1/4","1/2","1","2","4","8","16"))+
  scale_size("Number of mutants", range=c(1,8))+
  coord_flip(clip = "off")

library(patchwork)
scatterPlot=main_plot+customtheme
legend=ggpubr::as_ggplot(ggpubr::get_legend(main_plot))

export <-  scatterPlot + 
  legend + 
  plot_layout(nrow = 1, ncol=2,
              widths =  unit(c(8, 1), c('cm', 'null')),
              heights = unit(c(8, 1), c('cm', 'null')))

ggsave(
  plot = export,
  filename = "../mostDeviantHistologies6.pdf",
  device = cairo_pdf,
  width = 24,
  height = 24,
  units = "cm"
)

