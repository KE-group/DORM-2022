setwd("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Data/adding Targ to genomewide")
library(data.table)
rm(list=ls()); gc()

##------- Setting up data

# source('https://gist.githubusercontent.com/dchakro/8b1e97ba6853563dd0bb5b7be2317692/raw/parallelRDS.R')
# # Unused --------------
# # gwDT <- readRDS.gz("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/genome_wide/1_AllsamplesMinInfo.RDS")
# # 
# # gwDT <- gwDT[!grepl("_ENST", Gene.name, fixed = T), ]
# # gwDT <- gwDT[!stringi::stri_detect_regex(str = Mutation.AA, pattern = "^p.\\?"), ]
# # gwDT <- gwDT[!grepl(pattern = '=', x = Mutation.AA, fixed = T), ]
# # gwDT[, Mutation.AA := stringi::stri_replace_first_fixed(str = Mutation.AA,
# #                                                     pattern = "p.",
# #                                                     replacement = "")]
# # gwDT[, Mutation.AA := stringi::stri_replace_all_regex(str = Mutation.AA,
# #                                                   pattern = "\\*",
# #                                                   replacement = "X")]
# 
# fullDT <- readRDS.gz("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/targeted_and_wgs/1_AllsamplesMinInfo.RDS")
# 
# fullDT <- fullDT[!grepl("_ENST", Gene.name, fixed = T), ]
# fullDT <- fullDT[!stringi::stri_detect_regex(str = Mutation.AA, pattern = "^p.\\?"), ]
# fullDT <- fullDT[!grepl(pattern = '=', x = Mutation.AA, fixed = T), ]
# fullDT[, Mutation.AA := stringi::stri_replace_first_fixed(str = Mutation.AA,
#                                                         pattern = "p.",
#                                                         replacement = "")]
# fullDT[, Mutation.AA := stringi::stri_replace_all_regex(str = Mutation.AA,
#                                                       pattern = "\\*",
#                                                       replacement = "X")]
# 
# PMID_info <- unique(na.exclude(fullDT[,.(Sample.name,PMID,Gene.name)]))
# PubMedIDs <- unique(PMID_info$PMID)
# # Calculating the scale of study i.e. number of genes assayed in the individual studies
# number_of_genes <- PMID_info[,uniqueN (.SD),by=PMID,.SDcols=c("PMID","Gene.name")]
# rm(PMID_info);gc()
# 
# # Calculate the number of samples in the study (with atleast 1 mutation)
# countSamplesForPMID <- function(PubMedID){
#     return(nrow(unique(x = fullDT[PMID == PubMedID, "Sample.name"])))
# }
# 
# 
# results <- parallel::mclapply(X = PubMedIDs,
#                               FUN = countSamplesForPMID,
#                               mc.cores = parallel::detectCores())
# sampleCount_in_study <-data.table(PubMedID = PubMedIDs, N = unlist(results,use.names = F))


# ##--- Writing Data
# con <- pipe("pigz -p4 > 20221007_Data.GW.Targ.gz", "wb")
# save(fullDT, sampleCount_in_study, number_of_genes, file = con)
# close(con);rm(con)


