# Pies mutations normalized to number of samples.

rm(list=ls())
setwd("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/Pies/RAW")

muts <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v92/Full_Database/20210602_4.COSMIC.all.coding.Mutatations.RDS")
muts$Mutation.AA <- stringi::stri_replace_first_fixed(str = muts$Mutation.AA,pattern = "p.",replacement = "")
muts$Mutation.AA <- stringi::stri_replace_all_regex(str = muts$Mutation.AA,pattern = "\\*",replacement = "X")
# muts$mutID <- paste(muts$Gene.name,muts$Mutation.AA,sep="=")

# output <- plyr::count(muts,"mutID")
Stats <- as.data.frame(data.table::as.data.table(muts[,c("Gene.name","Sample.name")])[,.N, .(Gene.name,Sample.name)])
colnames(Stats)[3] <- "count"

Sample_Tissue_Map <- unique.data.frame(muts[,c("Sample.name","Primary.site")])
Stats$tissue <- Sample_Tissue_Map$Primary.site[match(x = Stats$Sample.name,table = Sample_Tissue_Map$Sample.name)]
Stats$tissue <- gsub("_"," ",Stats$tissue,fixed = T)

output <- as.data.frame(data.table::as.data.table(Stats[,c("Gene.name","tissue")])[,.N, .(Gene.name,tissue)])


rm(muts,Sample_Tissue_Map);gc()

library(ggplot2)

pie_theme <- theme_void() + theme(legend.position="right",
                                    legend.text=element_text(family = "serif",
                                                             size = 10),
                                    plot.title = element_text(family = "serif",
                                                              face="bold.italic",
                                                              color = "black"))

for(cType in sort(unique(tmp$Primary.site))){
  # Make a frequency table from the geneList
  sub <- data.frame(table(tmp[tmp$Primary.site == cType,2]))
  sub <- sub[order(sub[,2],decreasing = T),]
  colnames(sub) <- c("Gene","count")
  
  threshold <- 20
  pie_table_all <- aggregate(sub$count~sub$Gene,FUN=sum)
  colnames(pie_table_all) <- c("Gene","count")
  pie_table_all <- pie_table_all[order(pie_table_all$count,decreasing = T),]
  threshold <- min(threshold,length(pie_table_all$Gene))
  pie_table <- pie_table_all[1:threshold,]
  pie_table[(threshold+1),] <- list("Others",sum(pie_table_all$count[(threshold+1):length(pie_table_all$count)]))
  rm(pie_table_all,sub);gc()
  sliceColors <- rep(NA,length(pie_table$Gene))
  idx <- which(pie_table$Gene == "Others")
  sliceColors[idx] <- "#c7c7c7"
  sliceColors[-idx] <- viridis::plasma(length(pie_table$Gene[-idx]),direction = 1)
  names(sliceColors) <- pie_table$Gene
  print(cType)
  
  ggPie_1 <- ggplot(pie_table,aes(x="",y=count,fill=reorder(Gene,-count)))+geom_bar(stat="identity", width=1, color="white") + scale_fill_manual(values= sliceColors)+ guides(fill = guide_legend(title = "Genes", title.position = "top",title.theme = element_text(family="serif", face = "italic", angle = 0),override.aes = list(size = 3)))+ coord_polar(theta = "y",direction = -1)+ggtitle(paste(cType))+pie_theme
  
  ggplot(pie_table[-21,],aes(x="",y=count,fill=reorder(Gene,-count)))+geom_bar(stat="identity", width=1, color="white") + theme_void() + scale_fill_manual(values= sliceColors) + theme(legend.position="right",legend.text=element_text(family="serif",size=10)) + guides(fill = guide_legend(title = "Genes", title.position = "top",title.theme = element_text(family="serif", face = "italic", angle = 0),override.aes = list(size = 3)))+ coord_polar(theta = "y",direction = -1)
  print(pie_2)
  
  ggsave(filename = paste0("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/Pies/RAW/pie_",cType,".pdf"),
         plot = ggPie,
         width=6,
         height=5)
}
