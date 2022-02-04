#------> table_plyr_data.table------
# generate a table of instances of occurances of a string.
# e.g. TPCN1=M844I 1

test.name <- "table_plyr_data.table"
library(data.table)

source("https://gist.githubusercontent.com/dchakro/8b1e97ba6853563dd0bb5b7be2317692/raw/parallelRDS.R")
muts <- readRDS.gz("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Data/COSMIC_v92_R_DT/allCodingMutations.RDS")

muts$Mutation.AA <- stringi::stri_replace_first_fixed(str = muts$Mutation.AA,pattern = "p.",replacement = "")
muts$Mutation.AA <- stringi::stri_replace_all_regex(str = muts$Mutation.AA,pattern = "\\*",replacement = "X")
MutationID <- as.data.table(paste(muts$Gene.name, muts$Mutation.AA, sep = "="))
colnames(MutationID) <- c("MutationID")

rm(loadRDS, readRDS.gz, writeRDS, saveRDS.gz, muts)
gc()

source("https://raw.githubusercontent.com/dchakro/shared_Rscripts/master/summarySE.R")

DF <- data.frame(expr="",N=NA,time=NA,sd=NA,se=NA,ci=NA,size=NA,stringsAsFactors = F)
DF <- DF[-1,]

bmark_out_tsv <- paste0("table/bmark","_",test.name,".tsv")
file.create(bmark_out_tsv, showWarnings = F)
for(i in c(10000,1000000)){
  working_set <- MutationID[1:i,]
  bmark <- microbenchmark("table" = {
    var1 <- as.data.frame(base::table(working_set[,1]))
  },"plyr" = {
    var2 <- plyr::count(working_set,"MutationID")
  }, "data.tablex1"={
    setDTthreads(1)
    var3 <- working_set[,.N,.(MutationID)]
  },
  "data.tablex2"={
    setDTthreads(2)
    var4 <- working_set[,.N,.(MutationID)]
  }, times = 10)
  # saveRDS(bmark,file = paste0("./bmark/bmark_",test.name,"_",i,".RDS"))
  results <-
    summarySE(
      bmark,
      measurevar = "time",
      groupvars = "expr",
      statistic = "mean"
    )
  results$size <- rep(i,length(results[,1]))
  DF <- rbind.data.frame(DF,results)
  
  bmark$size <- rep(i,nrow(bmark))
  write.table(
    bmark,
    file = bmark_out_tsv,
    sep = "\t",
    row.names = F,
    col.names = F,
    append = T
  )
  rm(results,bmark)
}
# saveRDS(DF,file = paste0("./results/results_",test.name,".RDS"))

rm(list=ls())
gc()