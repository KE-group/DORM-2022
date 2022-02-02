library(ggplot2)
library(data.table)
rm(list=ls());gc()
# ---> Figure S2 <----

# Generating dataDF 
# source("https://gist.githubusercontent.com/dchakro/8b1e97ba6853563dd0bb5b7be2317692/raw/parallelRDS.R")
# muts <- readRDS.gz("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/Full_Database/4.COSMIC.all.coding.Mutatations.RDS")
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
# output <- muts[,.(.N), by = mutID]
# data.table::setnames(output, "N", "freq")
# singleOccurances <- output[freq==1, .(mutID)]
# rm(output)
# 
# '%nin%' <- Negate('%in%')
# dataDF <- muts[, .(Sample.name, mutID, Primary.site)]
# dataDF <- dataDF[mutID %nin% singleOccurances$mutID, ]
# rm("%nin%", singleOccurances, muts)
# gc()

# saveRDS.gz(dataDF, file="/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Data/COSMIC_v95_R_DT/RecurrentMuts.dataDF.RDS")

# Pre-load data
source("https://gist.githubusercontent.com/dchakro/8b1e97ba6853563dd0bb5b7be2317692/raw/parallelRDS.R")

dataDF <- readRDS.gz("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Data/COSMIC_v95_R_DT/RecurrentMuts.dataDF.RDS")
fullData <- readRDS.gz("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/Full_Database/1_AllsamplesMinInfo.RDS")

# <---- Purging irrelevant entries from data ---->
# Removing ENSEMBL transcripts and retaining just the genes
fullData <- fullData[!grepl("_ENST", Gene.name, fixed = T), ]
# Removing mutations annotated as ? in amino acid change
fullData <- fullData[!stringi::stri_detect_regex(str = Mutation.AA, pattern = "^p.\\?" ),]
# Removing Synonumous mutations
fullData <- fullData[!grepl(pattern = '=', x = Mutation.AA, fixed = T), ]

### Begin plotting

# ---> Figure S2 A : No. of samples with recurrent mutations by tissue type <----
rm(list = ls()[! ls() %in% c("dataDF","fullData")]) ; gc()
SampleDF <- unique(dataDF[,.(Sample.name, Primary.site)])
Stats <- SampleDF[,.(.N), by = Primary.site]
colnames(Stats) <- c("tissue","count")
setorder(Stats, -count)

SampleDF <- unique(fullData[,.(Sample.name, Primary.site)])
allStats <- SampleDF[,.(.N), by = Primary.site]
colnames(allStats) <- c("tissue","count")

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
ggplot(data = Stats, aes(x = reorder(tissue,-count),y = count,
                         fill=Type))+
  geom_col(width=0.75, position = "stack")+
  xlab("Tissue of origin of cancer")+
  ylab("Number of samples")+
  customtheme +
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = c("#c7c7c7","#000000"),
                    name="Type of \nSample")+
  theme(legend.position="right",
        legend.text=element_text(family="serif",
                                 size=11),
        legend.key.size = unit(0.5, "lines"))

ggsave(
  "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Figures/panels from R/S2 A.pdf",
  width = 8,
  height = 5,
  device = cairo_pdf
)


# ---> Figure S2 B : No. of Recurrent mutations by cancer type <----
rm(list = ls()[!ls() %in% c("dataDF", "fullData", "DC_theme_generator")])
gc()

Stats <- dataDF[,.(.N), by = Primary.site]
colnames(Stats) <- c("tissue","count")
Stats[, tissue := gsub("_", " ", tissue)]

allStats <- fullData[,.(.N), by = Primary.site]
colnames(allStats) <- c("tissue","count")

Stats <- merge.data.table(x = allStats, y = Stats, by = "tissue")
Stats[, tissue := gsub("_", " ", tissue)]
colnames(Stats)[2:3] <- c("All","Recurrent")
setorder(Stats,-All)

Stats[,Unique:=All-Recurrent]
summary(Stats$Recurrent/Stats$All)

Stats.bak <- Stats
Stats <- reshape2::melt(Stats[,c("tissue","Recurrent","Unique")])
colnames(Stats)[2:3] <- c("Type","count")
Stats$Type <- factor(Stats$Type,levels=c("Unique","Recurrent"))

if(!any(grepl("DC_theme_generator",x = ls()))){
  source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')  
}
customtheme <- DC_theme_generator(type = 'L',
                                  legend = 'F',
                                  ticks = 'out',
                                  x.axis.angle = 45,
                                  hjust = 1,
                                  vjust = 1,
                                  fontsize.cex = 1.2,
                                  ax.fontstyle = "italic")
options(scipen=100000)

ggplot(data = Stats, aes(x = reorder(tissue,-count),y = count,
                         fill=Type))+
  geom_col(width=0.75, position = "stack")+
  xlab("Tissue of origin of cancer")+
  ylab("Number of mutations")+
  customtheme +
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = c("#c7c7c7","#000000"),
                    name="Type of \nMutation")+
  theme(legend.position="right",
        legend.text=element_text(family="serif",
                                 size=11),
        legend.key.size = unit(0.5, "lines"))

ggsave(
  "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Figures/panels from R/S2 B.pdf",
  width = 8,
  height = 5,
  device = cairo_pdf
)

ggplot(data = Stats, aes(x = reorder(tissue,-count),y = count,
                         fill=Type))+
  geom_col(width=0.75, position = "fill")+
  geom_hline(yintercept = median(Stats.bak$Recurrent / Stats.bak$All), 
             color = "red",
             linetype = "dashed")+
  xlab("Tissue of origin of cancer")+
  ylab("Percentage of mutations")+
  customtheme +
  scale_y_continuous(expand = c(0, 0),
                     labels = paste0(c(0,25,50,75,100), " %"))+
  scale_fill_manual(values = c("#c7c7c7","#000000"),
                    name="Type of \nMutation")+
  theme(legend.position="right",
        legend.text=element_text(family="serif",
                                 size=11),
        legend.key.size = unit(0.5, "lines"))

ggsave(
  "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Figures/panels from R/S2 B_2.pdf",
  width = 8,
  height = 5,
  device = cairo_pdf
)


# ---> Figure S2 C : TMB plot <----
rm(list = ls()[! ls() %in% c("dataDF", "fullData", "DC_theme_generator")])

Stats <- dataDF[,.(.N), by = Sample.name]
colnames(Stats) <- c("Sample.name","count")
setorder(Stats, -count)
Sample_Tissue_Map <- unique(dataDF[,.(Sample.name, Primary.site)])
# print("Number of samples by tissue")
# Sample_Tissue_Map[, .N, .(Primary.site)]
Stats[, tissue := Sample_Tissue_Map$Primary.site[match(x = Stats$Sample.name, table = Sample_Tissue_Map$Sample.name)]]
Stats[, tissue := gsub("_", " ", tissue)]

if(!any(grepl("DC_theme_generator",x = ls()))){
  source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')  
}

customtheme <- DC_theme_generator(type = 'L',
                                  legend = 'F',
                                  ticks = 'out',
                                  x.axis.angle = 45,
                                  hjust = 1,
                                  vjust = 1,
                                  fontsize.cex = 1.2,
                                  ax.fontstyle = "italic")
options(scipen=100000)

ggplot(data = Stats, aes(y=count,
                         x=tissue))+
  geom_jitter(alpha=0.25,
              height = 0,
              size=1.5,
              width=0.3)+
  scale_y_continuous()+
  xlab("Tissue of origin of cancer")+
  ylab("Number of mutations\nper sample")+
  customtheme 

# ggsave(
#   "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/Barplots/TMB.pdf",
#   width = 8,
#   height = 5,
#   device = cairo_pdf
# )

ggsave(
  "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Figures/panels from R/S2 C.pdf",
  width = 8,
  height = 5,
  device = cairo_pdf
)
