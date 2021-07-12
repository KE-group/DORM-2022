rm(list=ls());gc()
library(data.table)
setwd("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/Normalized_to_sample")

plot_pie <- function(Tissue){
  # Tissue <- "lung"
  if(Tissue == "gastrointestinal tract (site indeterminate)"){
    Tissue_title <- "GI tract (site indeterminate)"
  }else if(Tissue == "haematopoietic and lymphoid tissue"){
    Tissue_title <- "haematopoietic / lymphoid"
  } else if(Tissue == "NS"){
    Tissue_title <- "Not specified"
  } else {
    Tissue_title <- Tissue
  }
  
  if(Tissue == "all"){
    subDF <- data.table::as.data.table(aggregate(DF$count ~ DF$Gene, FUN = sum))
    colnames(subDF) <- c("Gene.name","count")
    setorder(subDF,-count)
    plot_title <- paste0(Tissue_title," (n = ",sum(sampleCount$count),")")
  } else {
    subDF <- DF[tissue == Tissue ,]
    setorder(subDF,-count)
    subDF[,tissue:=NULL]
    plot_title <- paste0(Tissue_title," (n = ",sampleCount[tissue==Tissue,count],")")
  }
  # if slice weight normalized to total number of mutations in tissue type
  # subDF[, per_mutCount:=count/sum(subDF$count)]
  
  threshold <- 25
  threshold <- min(threshold, uniqueN(subDF))
  subDF <- subDF[1:threshold, ]
  '%nin%'=Negate('%in%')
  
  if(Tissue == "all"){
    tmp <- Stats 
  } else {
      tmp <- Stats[tissue == Tissue,]
  }
  
  samples_accounted_for <- tmp[Gene.name %in% unique(subDF$Gene.name),.(Sample.name)]
  samples_accounted_for <- unique(samples_accounted_for$Sample.name)
  others <- uniqueN(tmp[!Sample.name %in% samples_accounted_for,.(Sample.name)])
  rm(tmp,samples_accounted_for);gc()
  
  subDF <- data.table::rbindlist(l = list(subDF,list("Others",others)))
  if(Tissue != "all"){
    subDF[,per_sample:=count/sampleCount[tissue==Tissue,count]]  
  } else {
    subDF[,per_sample:=count/sum(sampleCount$count)]
  }
  
  pie_table <- subDF[,.(Gene.name,per_sample)]
  setnames(pie_table,"Gene.name","Gene")
  
  sliceColors <- rep(NA, uniqueN(pie_table))
  idx <- which(pie_table$Gene == "Others")
  sliceColors[idx] <- "#c7c7c7"
  sliceColors[-idx] <- viridis::plasma(uniqueN(pie_table)-1, direction = 1)
  names(sliceColors) <- pie_table$Gene
  pie_table$Gene <- factor(pie_table$Gene,levels = pie_table$Gene)
  print(paste(Tissue,round(sum(pie_table$per_sample),digits = 2)))
  
  library(ggplot2)
  pie <- ggplot(pie_table,aes(x="",
                              y=per_sample,
                              fill=Gene))+
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
          plot.title = element_text(hjust = 0.5),
          legend.key.size = unit(0.2, "lines")) + 
    guides(fill = guide_legend(title = "Genes", 
                               title.position = "top", 
                               byrow = T, 
                               nrow = 26, 
                               title.theme = element_text(family="serif", 
                                                          size = 6, 
                                                          face = "italic", 
                                                          angle = 0))) + 
    coord_polar(theta = "y",
                direction = -1) + 
    ggtitle(plot_title)
  
  return(pie)
}

# plot_pie("lung")
# plot_bar("pancreas")
# plot_bar("thyroid")

# ----> Set up environment <-------
source("https://gist.githubusercontent.com/dchakro/8b1e97ba6853563dd0bb5b7be2317692/raw/parallelRDS.R")

Stats <- readRDS.gz("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/COSMIC_v94_R_DT/CountStatsRAW.RDS")
rm(loadRDS,readRDS.gz,writeRDS,saveRDS.gz)

sampleCount <- unique(Stats[,.(Sample.name,tissue)])[,.N,.(tissue)]
setnames(sampleCount,"N","count")

DF <- Stats[,.N, .(Gene.name,tissue)]
setnames(DF,c("N"),c("count"))

tissues <- unique(DF$tissue)
# plot_pie("all")

myplots <-
  parallel::mclapply(
    X = as.list(c("all",sort(tissues))),
    FUN = function(X) plot_pie(Tissue = X),
    mc.cores = parallel::detectCores()
  )

ggplot2::ggsave(
  filename = "Pies_by_sample.pdf",
  plot = gridExtra::marrangeGrob(myplots, 
                                 layout_matrix = matrix(
                                   data = 1:16,
                                   nrow = 4,
                                   ncol = 4,
                                   byrow = T),
                                 as.table = F),
  width = 12,
  height = 14
)
