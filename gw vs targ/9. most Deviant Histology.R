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


library(ggplot2)
source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type = "L",legend = F)

main_plot <- ggplot(myDT,
       aes(
         x = FC,
         y = DIFF,
         color = tissue,
         size = full_mutant_N,
         alpha = gw_mutant_N
       )) +
  geom_point() +
  xlab("Fold change (+/-)\n(Full vs Genome-wide)")+
  ylab("Difference in Population Freq (%)\n(Full - Genome-wide)")+
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  scale_color_viridis_d(option = "turbo","Primary-Site") +
  scale_alpha("N (gw) / opacity",range = c(0.15, 0.8))+
  scale_size(range=c(1,8),trans = "log10")+
  coord_cartesian(clip = "off")

library(patchwork)
scatterPlot=main_plot+customtheme
legend=ggpubr::as_ggplot(ggpubr::get_legend(main_plot))

export <-  scatterPlot + 
  legend + 
  plot_layout(nrow = 1, ncol=2,
              widths =  unit(c(12, 1), c('cm', 'null')),
              heights = unit(c(12, 1), c('cm', 'null')))

ggsave(
  plot = export,
  filename = "../mostDeviantHistologies.pdf",
  device = cairo_pdf,
  width = 24,
  height = 24,
  units = "cm"
)

# plotDT <- myDT[(FC >=5 & FC < Inf)|(FC <= -5 & FC > -Inf),]
