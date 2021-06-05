rm(list=ls())

# muts <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v92/Full_Database/20210602_4.COSMIC.all.coding.Mutatations.RDS")
# muts$Mutation.AA <- stringi::stri_replace_first_fixed(str = muts$Mutation.AA,pattern = "p.",replacement = "")
# muts$Mutation.AA <- stringi::stri_replace_all_regex(str = muts$Mutation.AA,pattern = "\\*",replacement = "X")
# muts$mutID <- paste(muts$Gene.name,muts$Mutation.AA,sep="=")
# 
# output <- plyr::count(muts,"mutID")
# singleOccurances <- output$mutID[which(output$freq==1)]
# 
# '%nin%' <- Negate('%in%')
# tmp <- muts[muts$mutID %nin% singleOccurances,c("Sample.name","Gene.name","Primary.site")]
# rm("%nin%",singleOccurances,muts,output);gc()
# 
# saveRDS(object = tmp,file = "/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v92/Full_Database/20210605.Sample_Gene_pairs.RDS")


tmp <- unique.data.frame(readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v92/Full_Database/20210605.Sample_Gene_pairs.RDS"))

library(ggplot2)

pie_theme <- theme_void() + theme(legend.position="right",
                                    legend.text=element_text(family="serif",size=8),
                                    legend.key.height= unit(0.5, 'cm'),
                                    legend.key.width= unit(0.5, 'cm'),
                                    plot.title = element_text(family = "serif",
                                                              face="bold.italic",
                                                              color = "black"))

for(cType in sort(unique(tmp$Primary.site))){
  # Make a frequency table from the geneList
  sub <- data.frame(table(tmp[tmp$Primary.site == cType,2]))
  sub <- sub[order(sub[,2],decreasing = T),]
  
  threshold <- 20
  pie_table <- sub[1:threshold,]
  colnames(pie_table) <- c("Gene","count")
  pie_table$Gene <- as.character(pie_table$Gene)
  print(cType)
  
  ggPie <- ggplot(pie_table,aes(x="",y=count,fill=reorder(Gene,-count)))+geom_bar(stat="identity", width=1, color="white") + scale_fill_viridis_d(option = "plasma",direction = 1) + guides(fill = guide_legend(title = "Genes", title.position = "top",title.theme = element_text(family="serif", face = "italic", angle = 0,size=8)))+ coord_polar(theta = "y",direction = -1)+ggtitle(paste(cType))+pie_theme
  
  ggsave(filename = paste0("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/Pies/RAW/pie_",cType,".pdf"),
         plot = ggPie,
         width=6,
         height=5)
}
