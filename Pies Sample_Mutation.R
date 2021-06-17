# Pies mutations normalized to number of samples.
rm(list=ls());gc()
setwd("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/Pies/RAW")

# ----> Set up data <-------
# source("https://gist.githubusercontent.com/dchakro/8b1e97ba6853563dd0bb5b7be2317692/raw/parallelRDS.R")
# 
# # # muts <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v92/Full_Database/20210602_4.COSMIC.all.coding.Mutatations.RDS")
# 
# library(data.table)
# # # mutsDT <- data.table::as.data.table(muts)
# # # saveRDS.gz(object = mutsDT, file = "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/COSMIC_v92_R_DT/allCodingMutations.RDS")
# 
# muts <- readRDS.gz(file = "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/COSMIC_v92_R_DT/allCodingMutations.RDS")
# rm(loadRDS,readRDS.gz,writeRDS,saveRDS.gz)
# 
# muts$Mutation.AA <-
#   stringi::stri_replace_first_fixed(str = muts$Mutation.AA,
#                                     pattern = "p.",
#                                     replacement = "")
# muts$Mutation.AA <-
#   stringi::stri_replace_all_regex(str = muts$Mutation.AA,
#                                   pattern = "\\*",
#                                   replacement = "X")
# muts$mutID <- paste(muts$Gene.name, muts$Mutation.AA, sep = "=")
# 
# # output <- plyr::count(muts,"mutID")
# 
# paste("Total samples:", uniqueN(muts$Sample.name))
# paste("Total mutations:", uniqueN(muts$mutID))
# paste("Avg. mutation/sample:", uniqueN(muts$mutID) / uniqueN(muts$Sample.name))
# 
# 
# Stats <- muts[,.N, .(Gene.name,Sample.name)]
# # colnames(Stats)[3] <- "count"
# setnames(Stats,"N","count")
# 
# Sample_Tissue_Map <- unique(muts[,c("Sample.name","Primary.site")])
# print("Number of samples by tissue")
# Sample_Tissue_Map[, .N, .(Primary.site)]
# 
# Stats$tissue <- Sample_Tissue_Map$Primary.site[match(x = Stats$Sample.name,table = Sample_Tissue_Map$Sample.name)]
# Stats$tissue <- gsub("_"," ",Stats$tissue,fixed = T)
# rm(Sample_Tissue_Map,muts);gc()
# 
# # sampleCount <- Stats[,.N, .(tissue)]
# sampleCount <- unique(Stats[,.(Sample.name,tissue)])[,.N,.(tissue)]
# setnames(sampleCount,"N","count")
# saveRDS(object = sampleCount,file = "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/COSMIC_v92_R_DT/sampleCountByCancerType.RDS")
# 
# DF <- Stats[,.N, .(Gene.name,tissue)]
# setnames(DF,c("N"),c("count"))
# saveRDS(object = DF,file = "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/COSMIC_v92_R_DT/mutCountPerGeneByCancerType.RDS")
# rm(Stats);gc()

# --------> Data set up and saved <-------

plot_pie <- function(Tissue){
  # Tissue <- "lung"
  DF <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/COSMIC_v92_R_DT/mutCountPerGeneByCancerType.RDS")
  sampleCount <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/COSMIC_v92_R_DT/sampleCountByCancerType.RDS")
  if(Tissue == "all"){
    print("all")
    subDF <- data.table::as.data.table(aggregate(DF$count ~ DF$Gene, FUN = sum))
    colnames(subDF) <- c("Gene.name","count")
    subDF[,per_sample:=count/sum(sampleCount$count)]
  } else {
    subDF <- DF[tissue == Tissue ,]
    subDF[,per_sample:=count/sampleCount[tissue==Tissue,count]]
  }
  # slice weight normalized to number of samples for the tissue type
  
  
  # slice weight normalized to total number of mutations in tissue type
  # subDF[, per_mutCount:=count/sum(subDF$count)]
  
  pie_table_all <- subDF[,.(Gene.name,per_sample)]
  setnames(pie_table_all,"Gene.name","Gene")
  setorder(pie_table_all,-per_sample)
  
  threshold <- 25
  threshold <- min(threshold, uniqueN(pie_table_all))
  pie_table <- pie_table_all[1:threshold, ]
  # pie_table <- rbindlist(l = list(pie_table,list("Others",
  #                                     sum(pie_table_all[(threshold + 1):uniqueN(pie_table_all), .(per_sample)]))))
  rm(pie_table_all);gc()
  
  # sliceColors <- rep(NA, uniqueN(pie_table))
  # idx <- which(pie_table$Gene == "Others")
  # sliceColors[idx] <- "#c7c7c7"
  # sliceColors[-idx] <-
  sliceColors <- viridis::plasma(uniqueN(pie_table), direction = 1)
  names(sliceColors) <- pie_table$Gene
  
  library(ggplot2)
  # pie <- ggplot(pie_table,aes(x="",
  #                             y=per_sample,
  #                             fill=reorder(Gene,-per_sample)))+
  #   geom_bar(stat="identity", 
  #            width=1,  
  #            color=NA) + 
  #   theme_void() + 
  #   scale_fill_manual(values = sliceColors) + 
  #   theme(legend.position="right",
  #         legend.text=element_text(family="serif",
  #                                  size = 6),
  #         title = element_text(family = "serif", 
  #                              size = 7, 
  #                              face = "bold.italic"),
  #         plot.title = element_text(hjust = 0.5),
  #         legend.key.size = unit(0.2, "lines")) + 
  #   guides(fill = guide_legend(title = "Genes", 
  #                              title.position = "top", 
  #                              byrow = T, 
  #                              nrow = 21, 
  #                              title.theme = element_text(family="serif", 
  #                                                         size = 6, 
  #                                                         face = "italic", 
  #                                                         angle = 0))) + 
  #   coord_polar(theta = "y",
  #               direction = -1)+
  #   ggtitle(Tissue)
  
  bar <- ggplot(pie_table,aes(x=reorder(Gene,-per_sample),
                                     y=per_sample,
                                     fill=reorder(Gene,-per_sample)))+
    geom_bar(stat="identity", 
             position ="dodge", 
             width=0.75,  
             color=NA) + 
    scale_fill_manual(values = sliceColors) + 
    scale_y_continuous(limits = c(0,1),
                       expand = c(0,0), 
                       labels = paste0(seq(0,100,by=25),"%")) +
    ylab("Percentage of samples") + 
    xlab("Genes") +
    theme(axis.line = element_line(colour = "black",
                                   size=0.5),
          panel.border = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_line(colour = "black"),
          legend.position="right",
          legend.text=element_text(family="serif",
                                   size = 6),
          title = element_text(family = "serif", 
                               size = 8, 
                               face = "bold.italic"),
          plot.title = element_text(hjust = 0),
          axis.text.x = element_text(family = "serif",
                                     face = "plain",
                                     angle = 90,
                                     hjust = 1,
                                     vjust = 1,
                                     size = 6,
                                     color = "black"),
          axis.title = element_text(family="serif", 
                                    size = 8, 
                                    face = "italic", 
                                    angle = 0),
          legend.key.size = unit(0.2, "lines")) + 
    guides(fill = guide_legend(title = "Genes", 
                               title.position = "top", 
                               byrow = T, 
                               nrow = 25, 
                               title.theme = element_text(family="serif", 
                                                          size = 6, 
                                                          face = "italic", 
                                                          angle = 0))) +
    ggtitle(Tissue)
  
  return(bar)
    # ggtitle(paste0(gsub("_", " ", Title, fixed = T),
    #                " (n = ",
    #                prettyNum(
    #                  sum(pie_table$count,na.rm = T),
    #                  big.mark = " ",
    #                  scientific = F
    #                ),
    #                ")"))
  
}

plot_pie("lung")
plot_pie("all")

# subDF_FREQ$cent2 <- subDF_FREQ$N/sum(subDF_FREQ$N)
subDF_FREQ[, per_mutCount:=N/sum(subDF_FREQ$N)] # DT syntax


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
