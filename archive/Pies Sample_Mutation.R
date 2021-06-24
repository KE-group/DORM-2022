# Pies mutations normalized to number of samples.
plot_bar <- function(Tissue){
  # Tissue <- "lung"
  print(Tissue)
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
    subDF[,per_sample:=count/sum(sampleCount$count)]
    plot_title <- paste0(Tissue_title," (n = ",sum(sampleCount$count),")")
  } else {
    subDF <- DF[tissue == Tissue ,]
    subDF[,per_sample:=count/sampleCount[tissue==Tissue,count]]
    plot_title <- paste0(Tissue_title," (n = ",sampleCount[tissue==Tissue,count],")")
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
  
  # "Others" = samples with mutations in genes other than these 25
  rm(pie_table_all);gc()
  
  sliceColors <- viridis::plasma(uniqueN(pie_table), direction = 1)
  names(sliceColors) <- pie_table$Gene
  
  library(ggplot2)
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
          legend.position="none",
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
          axis.text.y = element_text(family = "serif",
                                     face = "plain",
                                     size = 6,
                                     color = "black"),
          axis.title = element_text(family="serif", 
                                    size = 8, 
                                    face = "italic", 
                                    angle = 0),
          legend.key.size = unit(0.2, "lines")) + 
    guides(fill = guide_legend(title = "Proteins", 
                               title.position = "top", 
                               byrow = T, 
                               nrow = 25, 
                               title.theme = element_text(family="serif", 
                                                          size = 6, 
                                                          face = "italic", 
                                                          angle = 0))) +
    ggtitle(plot_title)
  
  return(bar)
}

plot_pie <- function(Tissue){
  # Tissue <- "lung"
  print(Tissue)
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
  pie_table <- rbindlist(l = list(pie_table,list("Others",
                                      sum(pie_table_all[(threshold + 1):uniqueN(pie_table_all), .(per_sample)]))))
  rm(pie_table_all);gc()
  
  sliceColors <- rep(NA, uniqueN(pie_table))
  idx <- which(pie_table$Gene == "Others")
  sliceColors[idx] <- "#c7c7c7"
  sliceColors[-idx] <- viridis::plasma(uniqueN(pie_table)-1, direction = 1)
  names(sliceColors) <- pie_table$Gene
  
  library(ggplot2)
  pie <- ggplot(pie_table,aes(x="",
                              y=per_sample,
                              fill=reorder(Gene,-per_sample)))+
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
                               byrow = F,
                               title.theme = element_text(family="serif",
                                                          size = 6,
                                                          face = "italic",
                                                          angle = 0))) +
    coord_polar(theta = "y",
                direction = -1)+
    ggtitle(Tissue)
  
  return(pie)
    # ggtitle(paste0(gsub("_", " ", Title, fixed = T),
    #                " (n = ",
    #                prettyNum(
    #                  sum(pie_table$count,na.rm = T),
    #                  big.mark = " ",
    #                  scientific = F
    #                ),
    #                ")"))
}

DF <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/COSMIC_v92_R_DT/mutCountPerGeneByCancerType.RDS")
sampleCount <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/COSMIC_v92_R_DT/sampleCountByCancerType.RDS")

# plot_pie("lung")
# plot_bar("lung")
# 
# plot_bar("pancreas")
# plot_pie("pancreas")
# 
# plot_pie("all")
# plot_bar("all")


tissues <- unique(DF$tissue)
# plot_bar(Tissue = "all")

myplots <-
  parallel::mclapply(
    X = as.list(c("all",sort(tissues))),
    FUN = function(X) plot_bar(Tissue = X),
    mc.cores = parallel::detectCores()
  )

ggplot2::ggsave(
  filename = "Plots_combined.pdf",
  plot = gridExtra::marrangeGrob(myplots, 
                                 layout_matrix = matrix(
                                   data = 1:16,
                                   nrow = 4,
                                   ncol = 4,
                                   byrow = T),
                                 as.table = F),
  width = 12,
  height = 11
)
