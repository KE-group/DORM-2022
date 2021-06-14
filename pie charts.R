# Pies of mutations
rm(list = ls())
setwd("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/Pies/RAW")

DF <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/Git/GitHub/websites/eleniuslabtools.utu.fi/shiny-apps/COSMIC/HotspotMutations/data/20210603.FrequencyMutations.byTissue.RDS")
colnames(DF)[3] <- "all"


plot_pie <- function(subDF, Title){
  print(Title)
  subDF <- subDF[order(subDF[, 2], decreasing = T), ]
  colnames(subDF)[2] <- "count"
  threshold <- 20
  
  # Trying top N
  pie_table_all <- aggregate(subDF$count ~ subDF$Gene, FUN = sum)
  colnames(pie_table_all) <- c("Gene", "count")
  pie_table_all <-
    pie_table_all[order(pie_table_all$count, decreasing = T), ]
  threshold <- min(threshold, length(pie_table_all$Gene))
  pie_table <- pie_table_all[1:threshold, ]
  
  pie_table[(threshold + 1), ] <-
    list("Others", sum(pie_table_all$count[(threshold + 1):length(pie_table_all$count)]))
  rm(pie_table_all, subDF)
  gc()
  sliceColors <- rep(NA, length(pie_table$Gene))
  idx <- which(pie_table$Gene == "Others")
  sliceColors[idx] <- "#c7c7c7"
  sliceColors[-idx] <-
    viridis::plasma(length(pie_table$Gene[-idx]), direction = 1)
  names(sliceColors) <- pie_table$Gene
  
  library(ggplot2)
  ggplot(pie_table,aes(x="",
                       y=count,
                       fill=reorder(Gene,-count)))+
    geom_bar(stat="identity", 
             width=1, 
             color=NA) + 
    theme_void() + 
    scale_fill_manual(values = sliceColors) + 
    theme(legend.position="right",
          legend.text=element_text(family="serif",
                                   size = 6),
          title = element_text(family = "serif", 
                               size = 7, 
                               face = "bold.italic"),
          legend.key.size = unit(0.2, "lines")) + 
    guides(fill = guide_legend(title = "Genes", 
                               title.position = "top", 
                               byrow = T, 
                               nrow = 21, 
                               title.theme = element_text(family="serif", 
                                                          size = 6, 
                                                          face = "italic", 
                                                          angle = 0))) + 
    coord_polar(theta = "y",
                direction = -1)+
    ggtitle(paste0(gsub("_", " ", Title, fixed = T),
                   " (n = ",
                   prettyNum(
                     sum(pie_table$count),
                     big.mark = " ",
                     scientific = F
                   ),
                   ")"))
  
}

# lapply(
# X = c(3, 5),
# FUN = function(i)
#   plot_pie(subDF = na.omit(DF[, c(1, i)]), Title = colnames(DF)[i])
# )

myplots <-
  parallel::mclapply(
    X = c(3, 5:ncol(DF)),
    FUN = function(i)
      plot_pie(subDF = na.omit(DF[, c(1, i)]), Title = colnames(DF)[i]),
    mc.cores = parallel::detectCores()
  )


ggplot2::ggsave(filename = paste0("Pies_combined_1.pdf"),
       plot = gridExtra::grid.arrange(grobs=myplots[1:16], 
                                      nrow=4, 
                                      ncol=4),
       width = 14,
       height = 12)


ggplot2::ggsave(filename = paste0("Pies_combined_2.pdf"),
       plot = gridExtra::grid.arrange(grobs=myplots[17:32], 
                                      nrow=4, 
                                      ncol=4),
       width = 14,
       height = 12)


ggplot2::ggsave(filename = paste0("Pies_combined_3.pdf"),
       plot = gridExtra::grid.arrange(grobs=myplots[33:39], 
                                      nrow=4, 
                                      ncol=4),
       width = 14,
       height = 12)


