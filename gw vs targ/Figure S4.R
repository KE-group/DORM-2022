library(ggplot2)
library(data.table)
rm(list=ls());gc()
# ---> Figure S2 <----

# Generating dataDF 

# source("https://gist.githubusercontent.com/dchakro/8b1e97ba6853563dd0bb5b7be2317692/raw/parallelRDS.R")
# muts <- readRDS.gz("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/targeted_and_wgs/4.COSMIC.all.coding.Mutatations.RDS")
# 
# muts[, Mutation.AA := stringi::stri_replace_first_fixed(str = Mutation.AA,
#                                                         pattern = "p.",
#                                                         replacement = "")]
# 
# muts[, Mutation.AA := stringi::stri_replace_all_regex(str = Mutation.AA,
#                                                       pattern = "\\*",
#                                                       replacement = "X")]
# 
# muts[, mutID := paste(Gene.name,Mutation.AA,sep="=")]
# 
# keep <- c("Gene.name","Mutation.AA","Sample.name","Primary.site")
# fullData <- muts[, keep, with=F]
# saveRDS.gz(fullData, file="/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/targeted_and_wgs/AllMuts_fullData.RDS")
# 
# output <- muts[,.(.N), by = mutID]
# data.table::setnames(output, "N", "freq")
# singleOccurances <- output[freq==1, .(mutID)]
# rm(output)
# 
# '%nin%' <- Negate('%in%')
# dataDF <- muts[, .(Sample.name, mutID, Primary.site)]
# dataDF <- dataDF[mutID %nin% singleOccurances$mutID, ]
# saveRDS.gz(dataDF, file="/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/targeted_and_wgs/RecurrentMuts.dataDF.RDS")
# 
# rm("%nin%", singleOccurances, muts)
# gc()


# Pre-load data
source("https://gist.githubusercontent.com/dchakro/8b1e97ba6853563dd0bb5b7be2317692/raw/parallelRDS.R")

dataDF <- readRDS.gz("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/targeted_and_wgs/RecurrentMuts.dataDF.RDS")
fullData <- readRDS.gz("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/targeted_and_wgs/AllMuts_fullData.RDS")


### Begin plotting

# ---> Figure S4 A : No. of samples with recurrent mutations by tissue type <----
rm(list = ls()[! ls() %in% c("dataDF","fullData")]) ; gc()
SampleDF <- unique(dataDF[,.(Sample.name, Primary.site)])
Stats <- SampleDF[,.(.N), by = Primary.site]
colnames(Stats) <- c("tissue","count")
setorder(Stats, -count)

SampleDF_f <- unique(fullData[,.(Sample.name, Primary.site)])
allStats <- SampleDF_f[,.(.N), by = Primary.site]
colnames(allStats) <- c("tissue","count")
setorder(allStats, -count)

# rm(loadRDS,readRDS.gz,writeRDS,saveRDS.gz)

Stats <- merge.data.table(x = allStats, y = Stats, by = "tissue")
Stats[, tissue := gsub("_", " ", tissue)]
colnames(Stats)[2:3] <- c("All","Recurrent")
setorder(Stats,-All)
Stats[,Unique:=All-Recurrent]

Stats <- reshape2::melt(Stats[,c("tissue","Recurrent","Unique")])
colnames(Stats)[2:3] <- c("Type","count")
# Stats$Type <- relevel(Stats$Type,"Unique") # Makes Unique first
Stats$Type <- factor(Stats$Type,levels=c("Unique","Recurrent"))

source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type = 'L',
                                  legend = 'T',
                                  ticks = 'out',
                                  x.axis.angle = 45,
                                  hjust = 1,
                                  vjust = 1,
                                  fontsize.cex = 1.2,
                                  ax.fontstyle = "italic")
options(scipen=100000)
Stats$tissue <- reorder(Stats$tissue,-Stats$count)
ggplot(data = Stats, aes(x = tissue,
                         y = count,
                         fill=Type))+
  geom_col(width=0.75, position = "stack")+
  xlab("Tissue of origin of cancer")+
  ylab("Number of samples")+
  customtheme +
  scale_y_continuous(expand = c(0, 0),limits = c(0,max(Stats$count,na.rm = T)))+
  scale_fill_manual(values = c("#c7c7c7","#000000"),
                    name="Sample\n Harboring")+
  theme(legend.position="right",
        legend.text=element_text(family="serif",
                                 size=11),
        legend.key.size = unit(0.5, "lines"))



ggsave(
  "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Figures/panels from R/S4 A_gw.pdf",
  width = 8,
  height = 5,
  device = cairo_pdf
)
