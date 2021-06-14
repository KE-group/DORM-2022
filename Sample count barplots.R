# ---> Hotspot Mutations : Tisue and Mutation count stats <----

rm(list=ls())
library(ggplot2)
muts <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v92/Full_Database/20210602_4.COSMIC.all.coding.Mutatations.RDS")
muts$Mutation.AA <-
  stringi::stri_replace_first_fixed(str = muts$Mutation.AA,
                                    pattern = "p.",
                                    replacement = "")
muts$Mutation.AA <-
  stringi::stri_replace_all_regex(str = muts$Mutation.AA,
                                  pattern = "\\*",
                                  replacement = "X")
muts$mutID <- paste(muts$Gene.name, muts$Mutation.AA, sep = "=")

dataDF <- muts[, c("Sample.name", "mutID", "Primary.site")]
rm(muts)

tmpdf <-
  unique.data.frame(dataDF[, c("Sample.name", "Primary.site")])

Stats <- data.frame(table(tmpdf$Primary.site))
colnames(Stats) <- c("tissue", "count")
Stats$tissue <- as.character(Stats$tissue)
Stats$tissue <- gsub("_", " ", Stats$tissue, fixed = T)
Stats <- Stats[order(Stats$count, decreasing = T),]

source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(
  type = 'L',
  legend = 'F',
  ticks = 'out',
  x.axis.angle = 90,
  hjust = 1,
  vjust = 0.5,
  fontsize.cex = 1.2,
  ax.fontstyle = "italic"
)
options(scipen=100000)

ggplot(data = Stats,aes(y=count,
                        x=reorder(tissue,-count)))+
  geom_col(fill="#000000",
           width=0.75)+
  customtheme+
  scale_y_continuous(expand = c(0,0))+
  xlab("Tissue of origin of cancer")+
  ylab("Number of samples")

# ggsave("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/Barplots/SampleCounts.svg",width = 8,height = 5)

rm(tmpdf)
output <- plyr::count(dataDF, "mutID")
singleOccurances <- output$mutID[which(output$freq == 1)]
rm(output)

'%nin%' <- Negate('%in%')
dataDF <- dataDF[dataDF$mutID %nin% singleOccurances, ]
rm("%nin%", singleOccurances)

tmpdf <- unique.data.frame(dataDF[, c("Sample.name", "Primary.site")])

Stats <- data.frame(table(tmpdf$Primary.site))
colnames(Stats) <- c("tissue", "count")
Stats$tissue <- as.character(Stats$tissue)
Stats$tissue <- gsub("_", " ", Stats$tissue, fixed = T)
Stats <- Stats[order(Stats$count, decreasing = T), ]

options(scipen=100000)

ggplot(data = TissuesStats,
       aes(y=count,
           x=reorder(tissue,-count)))+
  geom_col(fill="#000000",width=0.75)+
  customtheme+
  scale_y_continuous(expand = c(0,0))+
  xlab("Tissue of origin of cancer")+
  ylab("Number of samples")

# ggsave("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/Barplots/SampleCounts_2.svg",width = 8,height = 5)
