library(data.table)
#-----------------------
# what does the targeted seq add to whole-genome data?
rm(list=ls()); gc()
source('https://gist.githubusercontent.com/dchakro/8b1e97ba6853563dd0bb5b7be2317692/raw/parallelRDS.R')
muts_gw <- readRDS.gz(file = "/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/Full_Database/4.COSMIC.all.coding.Mutatations.RDS",threads = 4)

muts_gw[, Mutation.AA := stringi::stri_replace_first_fixed(str = Mutation.AA,
                                                        pattern = "p.",
                                                        replacement = "")]

muts_gw[, Mutation.AA := stringi::stri_replace_all_regex(str = Mutation.AA,
                                                      pattern = "\\*",
                                                      replacement = "X")]

muts_gw[, mutID := paste(Gene.name,Mutation.AA,sep="=")]

SampleDF_gw <- unique(muts_gw[,.(Sample.name, Primary.site)]) ; nrow(SampleDF_gw) # 36 224

tmpdf <- muts_gw[,.(mutID,Primary.site)]
output <- tmpdf[,.N,.(mutID)]
data.table::setnames(output, "N", "freq")
singleOccurances_gw <- output[freq == 1, .(mutID)]; nrow(singleOccurances_gw) # 2 935 352
'%nin%' <- Negate('%in%')
tmpdf_gw <- tmpdf[mutID %nin% singleOccurances_gw$mutID,]; nrow(tmpdf_gw) # 1 887 757
rm("%nin%",tmpdf,output,singleOccurances_gw);gc()


muts_full <- readRDS.gz(file = "/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/targeted_and_wgs/4.COSMIC.all.coding.Mutatations.RDS",threads = 4)

muts_full[, Mutation.AA := stringi::stri_replace_first_fixed(str = Mutation.AA,
                                                           pattern = "p.",
                                                           replacement = "")]

muts_full[, Mutation.AA := stringi::stri_replace_all_regex(str = Mutation.AA,
                                                         pattern = "\\*",
                                                         replacement = "X")]

muts_full[, mutID := paste(Gene.name,Mutation.AA,sep="=")]

SampleDF_full <- unique(muts_full[,.(Sample.name, Primary.site)]) ; nrow(SampleDF_full) # 364 241

tmpdf <- muts_full[,.(mutID,Primary.site)]
output <- tmpdf[,.N,.(mutID)]
data.table::setnames(output, "N", "freq")
singleOccurances_full <- output[freq == 1, .(mutID)]; nrow(singleOccurances_full) # 3 029 924
'%nin%' <- Negate('%in%')
tmpdf_full <- tmpdf[mutID %nin% singleOccurances_full$mutID,]; nrow(tmpdf_full) # 2 369 282
rm("%nin%",tmpdf,output,singleOccurances_full);gc()


df_gw <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/Full_Database/20220117.FrequencyByMutation.RDS")
df_gw <- df_gw[,1:4]
df_gw[, MutID:=paste0(Gene,"=",Mutation)]

df_full <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/targeted_and_wgs/20220606.FrequencyByMutation.RDS")
df_full <- df_full[,1:4]
df_full[, MutID:=paste0(Gene,"=",Mutation)]

## Stats
# full has 12 % more muts than gw data
scales::label_percent()(abs((1-(nrow(muts_full)/nrow(muts_gw)))))

# full has 26 % more recurrent muts than gw data
scales::label_percent()(abs(1-(nrow(tmpdf_full)/nrow(tmpdf_gw))))

# full has 906 % more samples than gw data
scales::label_percent()(abs(1-(nrow(SampleDF_full)/nrow(SampleDF_gw))))

# full has 1184 % more mutations in top_100 than gw data
# View(df_full[which(df_full$MutID %in% df_gw$MutID[1:100]),])
N=100
C_gw <- sum(df_gw$counts[1:N]); C_gw # 17 169
C_full <- sum(df_full$counts[which(df_full$MutID %in% df_gw$MutID[1:100])]); C_full # 220 526
scales::label_percent()(abs(1-(C_full/C_gw)))

## Graph--


top100_full <- df_full[1:100,]
top100_gw <- df_gw[1:100,]

# intersect(top100_full$MutID,top100_gw$MutID)
# selection <- data.frame (MutID=unique(c(top100_full$MutID,top100_gw$MutID)))

selection <- data.frame(MutID=top100_gw$MutID)
selection$gw <- df_gw$counts[match(x = selection$MutID, table = df_gw$MutID)]
selection$full <- df_full$counts[match(x = selection$MutID, table = df_full$MutID)]
selection$full <- selection$full-ifelse(is.na(selection$gw),0,selection$gw)
selection$total <- selection$full+ifelse(is.na(selection$gw),0,selection$gw)

selection$gw_p <- (selection$gw/selection$total)*100
selection$full_p <- (selection$full/selection$total)*100


selection$MutID <- gsub("=", " ", selection$MutID)
write.table(
  x = selection,
  file = "../Data/adding Targ to genomewide/barplot.tsv",
  sep = "\t",
  row.names = F,
  col.names = T,
  quote = F
)
