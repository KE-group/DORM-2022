# --------> parallel saveRDS.gz() <-------
source("https://raw.githubusercontent.com/dchakro/shared_Rscripts/master/summarySE.R")
source('https://gist.githubusercontent.com/dchakro/8b1e97ba6853563dd0bb5b7be2317692/raw/parallelRDS.R')
DT <- readRDS.gz("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Data/COSMIC_v92_R_DT/allCodingMutations.RDS")

test.name <- "parallel_saveRDS"
size <- c("100", "10000", "1000000")

DF <- data.frame(expr="",N=NA,time=NA,sd=NA,se=NA,ci=NA,size=NA,stringsAsFactors = F)
DF <- DF[-1,] 

bmark_out_tsv <- paste0("table/bmark","_",test.name,".tsv")
file.create(bmark_out_tsv, showWarnings = F)
for(i in size){
  bmark <- microbenchmark(
    "base" = {
      base::saveRDS(object = DT[1:i],
                    compress = "gzip",
                    file = "/dev/null")
    },
    "parallel" = {
      saveRDS.gz(
        object = DT[1:i],
        threads = parallel::detectCores(),
        compression_level = 6,
        file = "/dev/null"
      )
    }
    ,times = 10)
  gc()
  # saveRDS(bmark,file = paste0("./bmark/bmark_",test.name,"_",i,".RDS"))
  results <- summarySE(bmark,measurevar = "time",groupvars = "expr",statistic = "mean")
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