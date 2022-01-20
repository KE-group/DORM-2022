rm(list=ls());gc()
setwd("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Data/COSMIC_v95_R_DT/")

# ----> Set up data <-------
source("https://gist.githubusercontent.com/dchakro/8b1e97ba6853563dd0bb5b7be2317692/raw/parallelRDS.R")

library(data.table)

# convert v92 to data.table 
# muts <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v92/Full_Database/20210602_4.COSMIC.all.coding.Mutatations.RDS")
# mutsDT <- data.table::as.data.table(muts)
# saveRDS.gz(object = mutsDT, file = "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/COSMIC_v92_R_DT/allCodingMutations.RDS")

muts <- readRDS.gz(file = "/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/Full_Database/4.COSMIC.all.coding.Mutatations.RDS")

muts$Mutation.AA <-
  stringi::stri_replace_first_fixed(str = muts$Mutation.AA,
                                    pattern = "p.",
                                    replacement = "")
muts$Mutation.AA <-
  stringi::stri_replace_all_regex(str = muts$Mutation.AA,
                                  pattern = "\\*",
                                  replacement = "X")
muts$mutID <- paste(muts$Gene.name, muts$Mutation.AA, sep = "=")

# output <- plyr::count(muts,"mutID")

paste("Total samples:", uniqueN(muts$Sample.name)) # 35 462
paste("Total mutations:", uniqueN(muts$mutID)) # 3 608 385
paste("Avg. mutation/sample:", uniqueN(muts$mutID) / uniqueN(muts$Sample.name)) # 101.75


Stats <- muts[,.N, .(Gene.name,Sample.name)]
# colnames(Stats)[3] <- "count"
setnames(Stats,"N","count")

Sample_Tissue_Map <- unique(muts[,.(Sample.name,Primary.site)])
print("Number of samples by tissue")
Sample_Tissue_Map[, .N, .(Primary.site)]

Stats$tissue <- Sample_Tissue_Map$Primary.site[match(x = Stats$Sample.name,table = Sample_Tissue_Map$Sample.name)]
Stats$tissue <- gsub("_"," ",Stats$tissue,fixed = T)
saveRDS.gz(object = Stats,file = "CountStatsRAW.RDS")
rm(Sample_Tissue_Map,muts);gc()

# sampleCount <- Stats[,.N, .(tissue)]
sampleCount <- unique(Stats[,.(Sample.name,tissue)])[,.N,.(tissue)]
setnames(sampleCount,"N","count")
saveRDS(object = sampleCount,file = "sampleCountByCancerType.RDS")

DF <- Stats[,.N, .(Gene.name,tissue)]
setnames(DF,c("N"),c("count"))
saveRDS(object = DF,file = "mutCountPerGeneByCancerType.RDS")
rm(loadRDS,readRDS.gz,writeRDS,saveRDS.gz)
rm(Stats);gc()

rm(list=ls());gc()
source("https://gist.githubusercontent.com/dchakro/8b1e97ba6853563dd0bb5b7be2317692/raw/parallelRDS.R")
dataDF <- readRDS.gz("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/Full_Database/1_AllsamplesMinInfo.RDS")

SampleDF <- unique(dataDF[,.(Sample.name, Primary.site)])
Stats <- SampleDF[,.(.N), by = Primary.site]
colnames(Stats) <- c("tissue","count")
setorder(Stats, -count)
saveRDS(object = Stats,file = "UnfilteredSampleCount.RDS")
rm(list=ls());gc()

# --------> Data set up and saved <-------
