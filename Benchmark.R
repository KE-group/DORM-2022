rm(list=ls())
# setwd("~/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/BenchmarkR/COSMIC_test/subset/")
setwd("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/Code optimization/")
dir.create("bmark")
dir.create("results")
library(microbenchmark)

#-----------------
## Reading TSV
# base, readr, vroom
source("https://raw.githubusercontent.com/dchakro/shared_Rscripts/master/summarySE.R")
test.name <- "readingTSV"
file.prefix <- c("101", "1001", "100001")
file.prefix <- paste0("~/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/BenchmarkR/COSMIC_test/subset/",file.prefix,"_d.tsv")
DF <-
  data.frame(
    expr = "",
    N = NA,
    time = NA,
    sd = NA,
    se = NA,
    ci = NA,
    size = NA,
    stringsAsFactors = F
  )
DF <- DF[-1,]
for(file.name in file.prefix){
  f <- tail(strsplit(file.name,"/")[[1]],n = 1)
  bmark <- microbenchmark(
    "base" = {
      dat <- utils::read.table(file = file.name,
                               header = T,
                               sep = "\t",
                               as.is = T,
                               stringsAsFactors = T)
      rm(dat)
    },
    "readr" = {
      dat <- readr::read_delim(file = file.name,
                               delim = "\t",
                               progress=F)
      rm(dat)
    },
    "vroom" = {
      dat <- vroom::vroom(file = file.name,
                          delim = "\t",
                          progress = F)
      rm(dat)
    }
    ,times = 10)
    gc()
    saveRDS(bmark,file = paste0("./bmark/bmark_",
                                test.name,
                                "_",
                                f,
                                ".RDS"))
    results <-
      summarySE(
        bmark,
        measurevar = "time",
        groupvars = "expr",
        statistic = "mean"
      )
    results$size <- rep(f, length(results[, 1]))
    DF <- rbind.data.frame(DF, results)
    rm(results, bmark)
}
saveRDS(DF,file = paste0("./results/results_",test.name,".RDS"))
rm(list=ls())
gc()

#----------------------------
# unlist & fixed
test.name <- "fixed_and_unlist"
dat <- utils::read.table(file = "~/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/BenchmarkR/COSMIC_test/subset/10001_d.tsv",header = T,sep = "\t",as.is = T,stringsAsFactors = T)

source("https://raw.githubusercontent.com/dchakro/shared_Rscripts/master/summarySE.R")

DF <- data.frame(expr="",N=NA,time=NA,sd=NA,se=NA,ci=NA,size=NA,stringsAsFactors = F)
DF <- DF[-1,]

field <- dat$Mutation.AA
for(i in c(100,1000,10000)){
  working_set <- field[1:i]
  bmark <- microbenchmark("default" = {
    res_1 <- base::substring(text = working_set,
                             first = {unlist(base::gregexpr(pattern = "p.",
                                                            text = working_set))+2})
  }, "faster"={
    res_2 <- base::substring(text = working_set,
                             first = {unlist(base::gregexpr(pattern = "p.",
                                                            text = working_set,
                                                            fixed = T),
                                             use.names = F)+2})
  }, times = 10,
  control = list("warmup"))
  # print(identical(res_1,res_2))
  saveRDS(bmark,file = paste0("./bmark/bmark_",test.name,"_",i,".RDS"))
  results <- summarySE(bmark,measurevar = "time",groupvars = "expr",statistic = "mean")
  results$size <- rep(i,length(results[,1]))
  DF <- rbind.data.frame(DF,results)
  rm(results,bmark)
}
saveRDS(DF,file = paste0("./results/results_",test.name,".RDS"))

rm(list=ls())
gc()

#----- data.frame vs data.table ------
test.name <- "data.frame_vs_data.table"
source("https://gist.githubusercontent.com/dchakro/8b1e97ba6853563dd0bb5b7be2317692/raw/parallelRDS.R")
library(data.table)

Stats <- readRDS.gz("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/COSMIC_v92_R_DT/CountStatsRAW.RDS")
rm(loadRDS,readRDS.gz,writeRDS,saveRDS.gz)

source("https://raw.githubusercontent.com/dchakro/shared_Rscripts/master/summarySE.R")

DF <- data.frame(expr="",N=NA,time=NA,sd=NA,se=NA,ci=NA,size=NA,stringsAsFactors = F)
DF <- DF[-1,]

for(i in c(100,10000,1000000)){
  working_set <- Stats[1:i]
  bmark <- microbenchmark("data.frame" = {
    var1 <- table(unique.data.frame(working_set[, c("Sample.name", "tissue")])[, "tissue"])
  }, "data.tablex1"={
    setDTthreads(1)
    var2 <- unique(working_set[,.(Sample.name,tissue)])[,.N,.(tissue)]
  },
  "data.tablex4"={
    setDTthreads(4)
    var3 <- unique(working_set[,.(Sample.name,tissue)])[,.N,.(tissue)]
  }, times = 10,
  control = list("warmup"))
  saveRDS(bmark,file = paste0("./bmark/bmark_",test.name,"_",i,".RDS"))
  results <- summarySE(bmark,measurevar = "time",groupvars = "expr",statistic = "mean")
  results$size <- rep(i,length(results[,1]))
  DF <- rbind.data.frame(DF,results)
  rm(results,bmark)
}
saveRDS(DF,file = paste0("./results/results_",test.name,".RDS"))

rm(list=ls())
gc()
#----------------------------
test.name <- "table_plyr_data.table"
library(data.table)

source("https://gist.githubusercontent.com/dchakro/8b1e97ba6853563dd0bb5b7be2317692/raw/parallelRDS.R")
muts <- readRDS.gz("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/COSMIC_v92_R_DT/allCodingMutations.RDS")

muts$Mutation.AA <- stringi::stri_replace_first_fixed(str = muts$Mutation.AA,pattern = "p.",replacement = "")
muts$Mutation.AA <- stringi::stri_replace_all_regex(str = muts$Mutation.AA,pattern = "\\*",replacement = "X")
Mutations <- paste(muts$Gene.name,muts$Mutation.AA,sep="=")
MutationID <- as.data.table(MutationID)

rm(loadRDS, readRDS.gz, writeRDS, saveRDS.gz, muts)
gc()

source("https://raw.githubusercontent.com/dchakro/shared_Rscripts/master/summarySE.R")

DF <- data.frame(expr="",N=NA,time=NA,sd=NA,se=NA,ci=NA,size=NA,stringsAsFactors = F)
DF <- DF[-1,]

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
  "data.tablex4"={
    setDTthreads(4)
    var4 <- working_set[,.N,.(MutationID)]
  }, times = 10,
  control = list("warmup"))
  saveRDS(bmark,file = paste0("./bmark/bmark_",test.name,"_",i,".RDS"))
  results <- summarySE(bmark,measurevar = "time",groupvars = "expr",statistic = "mean")
  results$size <- rep(i,length(results[,1]))
  DF <- rbind.data.frame(DF,results)
  rm(results,bmark)
}
saveRDS(DF,file = paste0("./results/results_",test.name,".RDS"))

rm(list=ls())
gc()
