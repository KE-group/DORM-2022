rm(list = ls())
DF <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/Git/GitHub/websites/eleniuslabtools.utu.fi/shiny-apps/COSMIC/HotspotMutations/DFa/20210603.FrequencyMutations.byTissue.RDS")

for(i in seq(5,ncol(DF))){
  sub <- na.omit(DF[,c(1,2,i)])
  sub <- sub[order(sub[,3],decreasing = T),]
  # Others cannot be represented as a slice
  threshold <- 20;sum(sub[1:threshold,3]);sum(sub[threshold+1:nrow(sub),3],na.rm = T)
  
  # Trying top N
  threshold <- 20
  topN <- sub[1:threshold,]
  colnames(topN) <- c("Gene","Mutation","count")
  pie_table <- aggregate(topN$count~topN$Gene,FUN=sum)
  colnames(pie_table) <- c("Gene","count")
  pie_table <- pie_table[order(pie_table[,2],decreasing = T),]
  
  # Trying to work with larger lists
  # pie_table[threshold+1,] <- list("others",sum(pie_table$count[10:nrow(pie_table)]))
  # pie_table <- pie_table[c(1:9,nrow(pie_table)),]
  
  library(ggplot2)
  pie <- ggplot(pie_table,aes(x="",y=count,fill=reorder(Gene,count)))+geom_bar(stat="identity", width=1, color="white") + theme_void() + scale_fill_viridis_d(option = "plasma",direction = -1) + theme(legend.position="bottom",legend.text=element_text(family="serif")) + guides(fill = guide_legend(title = "Genes", title.position = "left",title.theme = element_text(family="serif", face = "italic", angle = 0)))
  
  print(pie + coord_polar(theta = "y")) # same as # pie + coord_polar(theta = "y", start=0,direction = 1)
  
}
# pie + coord_polar()





rm(pie_table)
