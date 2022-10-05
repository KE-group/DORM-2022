rm(list=ls())
# setwd("~/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/BenchmarkR/COSMIC_test/subset/")
setwd("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Data/Code optimization/")
dir.create("bmark")
dir.create("results")
library(microbenchmark)

#------> fixed & unlist-----------
test.name <- "fixed_and_unlist"
dat <- data.table::fread(file = "~/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/BenchmarkR/COSMIC_test/subset/100001_d.tsv",header = T,sep = "\t")
colnames(dat) <- gsub(" ",".",colnames(dat))
source("https://raw.githubusercontent.com/dchakro/shared_Rscripts/master/summarySE.R")

DF <- data.frame(expr="",N=NA,time=NA,sd=NA,se=NA,ci=NA,size=NA,stringsAsFactors = F)
DF <- DF[-1,]

field <- dat$Mutation.AA
rm(dat); gc()
for(i in c(1000,100000)){
  working_set <- field[1:i]
  bmark <- microbenchmark("default" = {
    res_1 <- base::substring(text = working_set,
                             first = {unlist(base::gregexpr(pattern = "p.",
                                                            text = working_set))+2})
  }, "withParams"={
    res_2 <- base::substring(text = working_set,
                             first = {unlist(base::gregexpr(pattern = "p.",
                                                            text = working_set,
                                                            fixed = T),
                                             use.names = F)+2})
  }, times = 10)
  saveRDS(bmark, file = paste0("./bmark/bmark_", test.name, "_", i, ".RDS"))
  results <-
    summarySE(
      bmark,
      measurevar = "time",
      groupvars = "expr",
      statistic = "mean"
    )
  results$size <- rep(i,length(results[,1]))
  DF <- rbind.data.frame(DF,results)
  rm(results,bmark)
}
saveRDS(DF, file = paste0("./results/results_", test.name, ".RDS"))
rm(list=ls())
gc()