rm(list=ls())
library(ggplot2)
muts <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v92/Full_Database/20210602_4.COSMIC.all.coding.Mutatations.RDS")
muts$Mutation.AA <- stringi::stri_replace_first_fixed(str = muts$Mutation.AA,pattern = "p.",replacement = "")
muts$Mutation.AA <- stringi::stri_replace_all_regex(str = muts$Mutation.AA,pattern = "\\*",replacement = "X")
muts$mutID <- paste(muts$Gene.name,muts$Mutation.AA,sep="=")

dataDF <- muts[,c("Sample.name","mutID","Primary.site")]
rm(muts)

output <- plyr::count(dataDF,"mutID")
singleOccurances <- output$mutID[which(output$freq==1)]
rm(output)

'%nin%' <- Negate('%in%')
dataDF <- dataDF[dataDF$mutID %nin% singleOccurances,]
rm("%nin%",singleOccurances)

Stats <- data.frame(table(dataDF$Sample.name))
colnames(Stats) <- c("Sample.name","count")

Sample_Tissue_Map <- unique.data.frame(dataDF[,c("Sample.name","Primary.site")])

Stats$tissue <- Sample_Tissue_Map$Primary.site[match(x = Stats$Sample.name,table = Sample_Tissue_Map$Sample.name)]

source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type = 'L',
                                  legend = 'F',
                                  ticks = 'out',
                                  x.axis.angle = 90,
                                  hjust = 1,
                                  vjust = 0.5,
                                  fontsize.cex = 1.2,
                                  ax.fontstyle = "italic")
options(scipen=100000)
ggplot(data = Stats,aes(y=count,x=tissue))+geom_jitter(alpha=0.25)+scale_y_continuous()+xlab("Tissue of origin of cancer")+ylab("Tumor Mutation Burden")+customtheme

# +stat_summary(fun=mean, geom="crossbar", color="red",alpha=0.5)

# ggsave("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/Barplots/SampleCounts.svg",width = 8,height = 5)



# ------- Benchmarking -------

# library(microbenchmark)
# source("https://raw.githubusercontent.com/dchakro/shared_Rscripts/master/summarySE.R")
# 
# size <- 500000; bmark <- microbenchmark(
#   "plyr" = {
#     var1 <- plyr::count(dataDF[1:size,],"mutID")
#   },
#   "table" = {
#     var2 <- data.frame(table(dataDF$mutID[1:size]),stringsAsFactors = F)
#   },
#   data.table = { # Fastest
#     var3 <- data.table::as.data.table(dataDF[1:size,])[,.(.N), by = mutID]
#   }
#   ,times = 10)
# results <- summarySE(bmark,measurevar = "time",groupvars = "expr",statistic = "mean")
# 
# 
# var3 <- as.data.frame(var3)
# colnames(var2) <- colnames(var1)
# colnames(var3) <- colnames(var1)
# var2$mutID <- as.character(var2$mutID)
# 
# var1 <- var1[order(var1$freq,var1$mutID,decreasing = T),]
# var2 <- var2[order(var2$freq,var2$mutID,decreasing = T),]
# var3 <- var3[order(var3$freq,var3$mutID,decreasing = T),]
# 
# all(var1$freq==var2$freq)
# all(var1$mutID==var2$mutID)
# all(var1$freq==var3$freq)
# all(var1$mutID==var3$mutID)
