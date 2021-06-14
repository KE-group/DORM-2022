# ---> Hotspot Mutations : Tisue and Mutation count stats <----

rm(list=ls())
library(ggplot2)
muts <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v92/Full_Database/20210602_4.COSMIC.all.coding.Mutatations.RDS")
muts$Mutation.AA <- stringi::stri_replace_first_fixed(str = muts$Mutation.AA,pattern = "p.",replacement = "")
muts$Mutation.AA <- stringi::stri_replace_all_regex(str = muts$Mutation.AA,pattern = "\\*",replacement = "X")
muts$mutID <- paste(muts$Gene.name,muts$Mutation.AA,sep="=")

tmpdf <- muts[,c("mutID","Primary.site")]

TissuesStats <- data.frame(table(tmpdf$Primary.site))
colnames(TissuesStats) <- c("tissue","count")
TissuesStats$tissue <- as.character(TissuesStats$tissue)
TissuesStats$tissue <- gsub("_"," ",TissuesStats$tissue,fixed = T)
TissuesStats <- TissuesStats[order(TissuesStats$count,decreasing = T),]

source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type = 'L',
                                  legend = 'F',
                                  ticks = 'out',
                                  x.axis.angle = 45,
                                  hjust = 1,
                                  vjust = 1,
                                  fontsize.cex = 1.2,
                                  ax.fontstyle = "italic")
options(scipen=100000)
ggplot(data = TissuesStats,aes(y=count,x=reorder(tissue,-count)))+geom_col(fill="#000000",width=0.75)+customtheme+scale_y_continuous(expand = c(0,0))+xlab("Tissue of origin of cancer")+ylab("Number of mutations")
# ggsave("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/Barplots/CodingMuts.svg",width = 8,height = 5)

output <- plyr::count(tmpdf,"mutID")
singleOccurances <- output$mutID[which(output$freq==1)]

'%nin%' <- Negate('%in%')
tmpdf <- tmpdf[tmpdf$mutID %nin% singleOccurances,]
rm("%nin%",singleOccurances,output)

TissuesStats <- data.frame(table(tmpdf$Primary.site))
colnames(TissuesStats) <- c("tissue","count")
TissuesStats$tissue <- as.character(TissuesStats$tissue)
TissuesStats$tissue <- gsub("_"," ",TissuesStats$tissue,fixed = T)
TissuesStats <- TissuesStats[order(TissuesStats$count,decreasing = T),]

options(scipen=100000)
ggplot(data = TissuesStats,aes(y=count,x=reorder(tissue,-count)))+geom_col(fill="#000000",width=0.75)+customtheme+scale_y_continuous(expand = c(0,0))+xlab("Tissue of origin of cancer")+ylab("Number of mutations")
# ggsave("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/Barplots/Mutations_recurrent.svg",width = 8,height = 5)

# ---> Hotspot Residues : Tisue and Mutation count stats <----

rm(list=ls())
library(ggplot2)
muts <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v92/Full_Database/20210602_4.COSMIC.all.coding.Mutatations.RDS")
muts$Mutation.AA <- stringi::stri_replace_first_fixed(str = muts$Mutation.AA,pattern = "p.",replacement = "")
muts$Mutation.AA <- stringi::stri_replace_all_regex(str = muts$Mutation.AA,pattern = "\\*",replacement = "X")
muts$residue <- stringi::stri_extract_first(str = muts$Mutation.AA, regex = "[ACDEFGHIKLMNPQRSTVWYX]?[0-9]+(_[ACDEFGHIKLMNPQRSTVWYX]?[0-9]+)?")
muts$mutID <- paste(muts$Gene.name,muts$residue,sep="_")

tmpdf <- muts[,c("mutID","Primary.site")]
output <- plyr::count(tmpdf,"mutID")
singleOccurances <- output$mutID[which(output$freq==1)]

'%nin%' <- Negate('%in%')
tmpdf <- tmpdf[tmpdf$mutID %nin% singleOccurances,]
rm("%nin%",singleOccurances,output)

TissuesStats <- data.frame(table(tmpdf$Primary.site))
colnames(TissuesStats) <- c("tissue","count")
TissuesStats$tissue <- as.character(TissuesStats$tissue)
TissuesStats$tissue <- gsub("_"," ",TissuesStats$tissue,fixed = T)
TissuesStats <- TissuesStats[order(TissuesStats$count,decreasing = T),]

source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type = 'L',
                                  legend = 'F',
                                  ticks = 'out',
                                  x.axis.angle = 45,
                                  hjust = 1,
                                  vjust = 1,
                                  fontsize.cex = 1.2,
                                  ax.fontstyle = "italic")

options(scipen=100000)
ggplot(data = TissuesStats,aes(y=count,x=reorder(tissue,-count)))+geom_col(fill="#000000",width=0.75)+customtheme+scale_y_continuous(expand = c(0,0))+xlab("Tissue of origin of cancer")+ylab("Number of mutations")
# ggsave("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/Barplots/Residues_recurrent.svg",width = 8,height = 5)
