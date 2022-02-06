#------> Reading TSV <--------
# base, readr, vroom
source("https://raw.githubusercontent.com/dchakro/shared_Rscripts/master/summarySE.R")
test.name <- "readingTSV"
# file.prefix <- c("101", "1001", "100001")
file.prefix <- c("1001", "100001")
file.prefix <- paste0("~/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/BenchmarkR/COSMIC_test/subset/",file.prefix,"_d.tsv")
bmark_out_tsv <- paste0("table/bmark","_",test.name,".tsv")
file.create(bmark_out_tsv, showWarnings = F)
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
    },
    "data.table" = {
      dat <- data.table::fread(file = file.name,
                               sep = "\t", 
                               header = T, 
                               showProgress = F)
      rm(dat)
    }
    ,times = 10)
  gc()
  # saveRDS(bmark, file = paste0("./bmark/bmark_", test.name, "_", f, ".RDS"))
  results <-
    summarySE(
      bmark,
      measurevar = "time",
      groupvars = "expr",
      statistic = "mean"
    )
  results$size <- rep(f, length(results[, 1]))
  DF <- rbind.data.frame(DF, results)
  bmark$size <- rep(f,nrow(bmark))
  write.table(
    bmark,
    file = bmark_out_tsv,
    sep = "\t",
    row.names = F,
    col.names = F,
    append = T
  )
  rm(results, bmark)
}
# saveRDS(DF,file = paste0("./results/results_",test.name,".RDS"))
rm(list=ls())
gc()