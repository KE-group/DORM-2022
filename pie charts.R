# Pies of mutations
rm(list = ls())
setwd("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/Pies/RAW")

DF <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/Git/GitHub/websites/eleniuslabtools.utu.fi/shiny-apps/COSMIC/HotspotMutations/data/20210603.FrequencyMutations.byTissue.RDS")
colnames(DF)[3] <- "all"


plot_pie <- function(subDF, Title){
  print(Title)
  subDF <- subDF[order(subDF[,2],decreasing = T),]
  colnames(subDF)[2] <- "count"
  threshold <- 20
  
  # Trying top N
  pie_table_all <- aggregate(subDF$count~subDF$Gene,FUN=sum)
  colnames(pie_table_all) <- c("Gene","count")
  pie_table_all <- pie_table_all[order(pie_table_all$count,decreasing = T),]
  threshold <- min(threshold,length(pie_table_all$Gene))
  pie_table <- pie_table_all[1:threshold,]

    pie_table[(threshold+1),] <- list("Others",sum(pie_table_all$count[(threshold+1):length(pie_table_all$count)]))
  rm(pie_table_all,subDF);gc()
  sliceColors <- rep(NA,length(pie_table$Gene))
  idx <- which(pie_table$Gene == "Others")
  sliceColors[idx] <- "#c7c7c7"
  sliceColors[-idx] <- viridis::plasma(length(pie_table$Gene[-idx]),direction = 1)
  names(sliceColors) <- pie_table$Gene
  
  library(ggplot2)
  ggplot(pie_table,aes(x="",y=count,fill=reorder(Gene,-count)))+geom_bar(stat="identity", width=1, color=NA) + theme_void() + scale_fill_manual(values= sliceColors) + theme(legend.position="right",legend.text=element_text(family="serif",size = 8),title = element_text(family = "serif", size = 14, face = "bold.italic"),legend.key.size = unit(0.2, "lines")) + guides(fill = guide_legend(title = "Genes", title.position = "top", byrow = T, nrow = (threshold+1), title.theme = element_text(family="serif", size = 10, face = "italic", angle = 0))) + coord_polar(theta = "y",direction = -1)+ggtitle(paste0(gsub("_"," ",Title,fixed=T)," (n = ",prettyNum(sum(pie_table$count),big.mark = " ",scientific = F),")"))
  
}

myplots <- lapply(X = c(3,5),FUN = function(i) plot_pie(subDF = na.omit(DF[,c(1,i)]), Title = colnames(DF)[i]))

myplots <- lapply(X = c(3,5:ncol(DF)),FUN = function(i) plot_pie(subDF = na.omit(DF[,c(1,i)]), Title = colnames(DF)[i]))

tmp <- gridExtra::grid.arrange(grobs=myplots, 
                               nrow=3, 
                               ncol=13)
ggsave(filename = paste0("Pies_combined.svg"),plot = tmp,width = 20,height = 10)







pdf("pies.pdf",height = 5,width = 6,onefile = T)

for(i in c(3,5:ncol(DF))){
  # i <- 21

  
  
  library(ggplot2)
   pie_1 <- ggplot(pie_table[-21,],aes(x="",y=count,fill=reorder(Gene,-count)))+geom_bar(stat="identity", width=1, color="white") + theme_void() + scale_fill_manual(values= sliceColors) + theme(legend.position="none") + coord_polar(theta = "y",direction = -1)
  
  
  
  pie <- gridExtra::grid.arrange(pie_1,
                           pie_2,
                           ncol=2,
                           nrow=1,
                           widths = c(2, 5),
                           top = grid::textGrob(label = paste0(gsub("_"," ",colnames(DF)[i],fixed=T)," (n = ",prettyNum(sum(pie_table$count),big.mark = " ",scientific = F),")"), gp = grid::gpar(fontsize=18, font=4, fontfamily="serif"),just = "right",y = -1))
  
  print(pie)
  
} ; dev.off()
