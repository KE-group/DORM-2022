rm(list=ls());gc()
setwd("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/Pies/Normalized_to_sample/")

# ----> Set up data <-------
source("https://gist.githubusercontent.com/dchakro/8b1e97ba6853563dd0bb5b7be2317692/raw/parallelRDS.R")

# muts <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v92/Full_Database/20210602_4.COSMIC.all.coding.Mutatations.RDS")

library(data.table)
setDTthreads(4)
# mutsDT <- data.table::as.data.table(muts)
# saveRDS.gz(object = mutsDT, file = "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/COSMIC_v92_R_DT/allCodingMutations.RDS")

muts <- readRDS.gz(file = "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/COSMIC_v92_R_DT/allCodingMutations.RDS")

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

paste("Total samples:", uniqueN(muts$Sample.name))
paste("Total mutations:", uniqueN(muts$mutID))
paste("Avg. mutation/sample:", uniqueN(muts$mutID) / uniqueN(muts$Sample.name))


Stats <- muts[,.N, .(Gene.name,Sample.name)]
# colnames(Stats)[3] <- "count"
setnames(Stats,"N","count")

Sample_Tissue_Map <- unique(muts[,c("Sample.name","Primary.site")])
print("Number of samples by tissue")
Sample_Tissue_Map[, .N, .(Primary.site)]

Stats$tissue <- Sample_Tissue_Map$Primary.site[match(x = Stats$Sample.name,table = Sample_Tissue_Map$Sample.name)]
Stats$tissue <- gsub("_"," ",Stats$tissue,fixed = T)
saveRDS.gz(object = Stats,file = "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/COSMIC_v92_R_DT/CountStatsRAW.RDS")
rm(Sample_Tissue_Map,muts);gc()

# sampleCount <- Stats[,.N, .(tissue)]
sampleCount <- unique(Stats[,.(Sample.name,tissue)])[,.N,.(tissue)]
setnames(sampleCount,"N","count")
saveRDS(object = sampleCount,file = "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/COSMIC_v92_R_DT/sampleCountByCancerType.RDS")

DF <- Stats[,.N, .(Gene.name,tissue)]
setnames(DF,c("N"),c("count"))
saveRDS(object = DF,file = "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/COSMIC_v92_R_DT/mutCountPerGeneByCancerType.RDS")
rm(loadRDS,readRDS.gz,writeRDS,saveRDS.gz)
rm(Stats);gc()


# --------> Data set up and saved <-------
