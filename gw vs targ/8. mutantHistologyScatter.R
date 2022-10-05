setwd("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Data/adding Targ to genomewide/mutants/")
library(data.table)
rm(list=ls()); gc()

#----------
files <- setdiff(list.files(), list.dirs(recursive = FALSE, full.names = FALSE))
        