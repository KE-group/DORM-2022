# -------> parallel readRDS.gz() <-------
source("https://raw.githubusercontent.com/dchakro/shared_Rscripts/master/summarySE.R")
source('https://gist.githubusercontent.com/dchakro/8b1e97ba6853563dd0bb5b7be2317692/raw/parallelRDS.R')
FILE <- "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Data/COSMIC_v95_R_DT/CountStatsRAW.RDS" # 51 MB RDS file

test.name <- "parallel_readRDS"

DF <- data.frame(expr="",N=NA,time=NA,sd=NA,se=NA,ci=NA,size=NA,stringsAsFactors = F)
DF <- DF[-1,] 

bmark_out_tsv <- paste0("table/bmark","_",test.name,".tsv")
file.create(bmark_out_tsv, showWarnings = F)
bmark <- microbenchmark(
  "base" = {
    base::readRDS(file = FILE)
  },
  "parallel" = {
    readRDS.gz(file = FILE, threads = parallel::detectCores())
  }
  ,times = 10)
gc()
# saveRDS(bmark,file = paste0("./bmark/bmark_",test.name,".RDS"))
results <- summarySE(bmark,measurevar = "time",groupvars = "expr",statistic = "mean")
results$size <- rep("51 MB", length(results[, 1]))
DF <- rbind.data.frame(DF, results)

bmark$size <- rep("51 MB",nrow(bmark))
write.table(
  bmark,
  file = bmark_out_tsv,
  sep = "\t",
  row.names = F,
  col.names = F,
  append = T
)
rm(results, bmark)
# saveRDS(DF, file = paste0("./results/results_", test.name, ".RDS"))
rm(list = ls())
gc()

