# Pies of all mutations (by tissue type)
rm(list = ls())
library(data.table)
setwd("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/Pies/RAW")

DF <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/Git/GitHub/websites/eleniuslabtools.utu.fi/shiny-apps/COSMIC/HotspotMutations/data/20210603.FrequencyMutations.byTissue.RDS")
colnames(DF)[3] <- "all"


plot_pie <- function(subDT){
  subDT <- data.table:::na.omit.data.table( data.table::as.data.table(subDT))
  Title <- colnames(subDT)[2]
  print(Title)
  setorderv(subDT, colnames(subDT)[2],order = -1)
  colnames(subDT)[2] <- "count"
  threshold <- 20
  
  # Trying top N
  pie_table_all <- subDT[,.(x=sum(count)), by = Gene]
  data.table::setnames(pie_table_all, old = "x", new = "count")
  data.table::setorder(pie_table_all, -count, Gene)
  # pie_table_all <- aggregate(subDT$count ~ subDT$Gene, FUN = sum)
  # colnames(pie_table_all) <- c("Gene", "count")
  # pie_table_all <-
  # pie_table_all[order(pie_table_all$count, decreasing = T), ]
  
  threshold <- min(threshold, length(pie_table_all$Gene))
  pie_table <- pie_table_all[1:threshold, ]
  
  others <-
    sum(pie_table_all[(threshold + 1):dim(pie_table_all)[1], "count"], na.rm = T)
  pie_table <-
    data.table::rbindlist(l = list(pie_table, list("Others", others)))
  
  # pie_table[(threshold + 1), ] <-
  #   list("Others", sum(pie_table_all$count[(threshold + 1):length(pie_table_all$count)]))
  rm(pie_table_all, subDT)
  gc()
  sliceColors <- rep(NA, length(pie_table$Gene))
  idx <- which(pie_table$Gene == "Others")
  sliceColors[idx] <- "#c7c7c7"
  sliceColors[-idx] <-
    viridis::plasma(length(pie_table$Gene[-idx]), direction = 1)
  names(sliceColors) <- pie_table$Gene
  
  library(ggplot2)
  pie <- ggplot(pie_table,aes(x="",
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
          plot.title = element_text(hjust = 0.5),
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
                     sum(pie_table$count,na.rm = T),
                     big.mark = " ",
                     scientific = F
                   ),
                   ")"))
  return(pie)
  
}

# lapply(X = c(3,21),
#        FUN = function(i) plot_pie(subDT = DF[, c(1, i)]))



# Alphabetical order
myplots <-
  parallel::mclapply(
    X = c(3, 5:ncol(DF)),
    FUN = function(i)
      plot_pie(subDT = DF[, c(1, i)]),
    mc.cores = parallel::detectCores()
  )

# Ordered by total number of mutations
# No_of_mutations <- colSums(DF[, c(5:ncol(DF))], na.rm = T)
# range(order(No_of_mutations,decreasing = T)+4)
# 
# myplots <-
#   parallel::mclapply(
#     X = c(3, (order(No_of_mutations, decreasing = T)+4)),
#     FUN = function(i)
#       plot_pie(subDT = DF[, c(1, i)]),
#     mc.cores = parallel::detectCores()
#   )

ggplot2::ggsave(
  filename = "Pies_combined.pdf",
  plot = gridExtra::marrangeGrob(myplots, 
                                 layout_matrix = matrix(
                                   data = 1:16,
                                   nrow = 4,
                                   ncol = 4,
                                   byrow = T),
                                 as.table = F),
  width = 14,
  height = 12
)

