rm(list=ls());gc()
# setwd("~/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/BenchmarkR/COSMIC_test/subset/")
setwd("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Data/Code optimization/")
library(microbenchmark)

dir.create("bmark")
dir.create("results")
dir.create("table")

#------> Reading TSV <--------
rm(list=ls());gc()
source("../../code/benchmark_snippet/reading_TSV.R")

#------> data.frame vs data.table ------
rm(list=ls());gc()
source("../../code/benchmark_snippet/DF vs DT.R")

#------> table_plyr_data.table------
rm(list=ls());gc()
source("../../code/benchmark_snippet/table_vs_plyr_vs_data.table.R")

#------> base vs stringi <-------
rm(list=ls());gc()
source("../../code/benchmark_snippet/base_vs_stringi.R")

#------> parallel processing <------
rm(list=ls());gc()
source("../../code/benchmark_snippet/parallel processing.R")

# --------> parallel saveRDS.gz() <-------
rm(list=ls());gc()
source("../../code/benchmark_snippet/parallel_saveRDS.R")

# -------> parallel readRDS.gz() <-------
rm(list=ls());gc()
source("../../code/benchmark_snippet/parallel_readRDS.R")

# --------> grep: GNU vs BSD vs R <--------
rm(list=ls());gc()
source("../../code/benchmark_snippet/grep.R")
