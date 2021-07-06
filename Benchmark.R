rm(list=ls())
# setwd("~/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/BenchmarkR/COSMIC_test/subset/")
setwd("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/Code optimization/")
dir.create("bmark")
dir.create("results")
library(microbenchmark)

#------> Reading TSV <--------
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
    rm(results, bmark)
}
# saveRDS(DF,file = paste0("./results/results_",test.name,".RDS"))
rm(list=ls())
gc()

#------> fixed & unlist-----------
test.name <- "fixed_and_unlist"
dat <- data.table::fread(file = "~/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/BenchmarkR/COSMIC_test/subset/10001_d.tsv",header = T,sep = "\t")
colnames(dat) <- gsub(" ",".",colnames(dat))
source("https://raw.githubusercontent.com/dchakro/shared_Rscripts/master/summarySE.R")

DF <- data.frame(expr="",N=NA,time=NA,sd=NA,se=NA,ci=NA,size=NA,stringsAsFactors = F)
DF <- DF[-1,]

field <- dat$Mutation.AA
rm(dat); gc()
for(i in c(100,1000,10000)){
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
  # print(identical(res_1,res_2))
  # saveRDS(bmark, file = paste0("./bmark/bmark_", test.name, "_", i, ".RDS"))
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
# saveRDS(DF, file = paste0("./results/results_", test.name, ".RDS"))

rm(list=ls())
gc()

#------> data.frame vs data.table ------
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
  }, times = 10)
  # saveRDS(bmark, file = paste0("./bmark/bmark_", test.name, "_", i, ".RDS"))
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
# saveRDS(DF, file = paste0("./results/results_", test.name, ".RDS"))

rm(list=ls())
gc()


#------> table_plyr_data.table------
test.name <- "table_plyr_data.table"
library(data.table)

source("https://gist.githubusercontent.com/dchakro/8b1e97ba6853563dd0bb5b7be2317692/raw/parallelRDS.R")
muts <- readRDS.gz("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/COSMIC_v92_R_DT/allCodingMutations.RDS")

muts$Mutation.AA <- stringi::stri_replace_first_fixed(str = muts$Mutation.AA,pattern = "p.",replacement = "")
muts$Mutation.AA <- stringi::stri_replace_all_regex(str = muts$Mutation.AA,pattern = "\\*",replacement = "X")
MutationID <- as.data.table(paste(muts$Gene.name, muts$Mutation.AA, sep = "="))
colnames(MutationID) <- c("MutationID")

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
  "data.tablex2"={
    setDTthreads(2)
    var4 <- working_set[,.N,.(MutationID)]
  }, times = 10)
  saveRDS(bmark,file = paste0("./bmark/bmark_",test.name,"_",i,".RDS"))
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
saveRDS(DF,file = paste0("./results/results_",test.name,".RDS"))

rm(list=ls())
gc()

#------> base vs stringi <-------
test.name <- "base_V_stringi"
dat <-
  data.table::fread(file = "~/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/BenchmarkR/COSMIC_test/subset/1000001_d.tsv", header = T, sep = "\t")
colnames(dat) <- gsub(" ",".",colnames(dat))

dat$Mutation.AA <- gsub("p.", "", dat$Mutation.AA, fixed = T)
dat <-
  dat[-grep("?", dat$Mutation.AA, fixed = T), ] # Removing "Unknown" mutations
dat <-
  dat[-grep("=", dat$Mutation.AA, fixed = T), ] # Removing silent mutations

source("https://raw.githubusercontent.com/dchakro/shared_Rscripts/master/summarySE.R")
DF <- data.frame(expr="",N=NA,time=NA,sd=NA,se=NA,ci=NA,size=NA,stringsAsFactors = F)
DF <- DF[-1,]
regex_pattern <- "[ACDEFGHIKLMNPQRSTVWYX]?[0-9]+"

for(i in c(100,10000,length(dat$Mutation.AA))){
  Mutation.AA <- dat$Mutation.AA[1:i]
  bmark <- microbenchmark("base" = {
    # res_b <- find_length_base(genomePos = genomePosition)
    res_b <-
      unlist(lapply(base::regmatches(
        x = Mutation.AA,
        m = gregexpr(pattern = regex_pattern, text = Mutation.AA, fixed = F)
      ), `[[`, 1))
    
  }, "stringi"={
    res_s <- stringi::stri_extract_first(str = Mutation.AA, regex = regex_pattern)
  }, times = 10)
  # print(identical(res_b,res_s))
  saveRDS(bmark,file = paste0("./bmark/bmark_",test.name,"_",i,".RDS"))
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
saveRDS(DF,file = paste0("./results/results_",test.name,".RDS"))
rm(list=ls())
gc()

# parallel processing
#-------->for vs apply vs foreach vs parSapply vs mapply <----------

test.name <- "serialVparallel"
library(doParallel)
library(data.table)
source("https://raw.githubusercontent.com/dchakro/shared_Rscripts/master/summarySE.R")
dat <-
  data.table::fread(file = "~/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/BenchmarkR/COSMIC_test/subset/1000001_d.tsv", header = T, sep = "\t")
colnames(dat) <- gsub(" ",".",colnames(dat))

dat <- dat[-grep("_ENST", dat$Gene.name, fixed = T), ]
dat$Mutation.AA <- gsub("p.", "", dat$Mutation.AA, fixed = T)
dat <- dat[-grep("?", dat$Mutation.AA, fixed = T), ]
dat <- dat[-grep("=", dat$Mutation.AA, fixed = T), ]
# dim(dat) # 106450     40

dat$mutID <- paste(dat$Gene.name, dat$Mutation.AA, sep = "=")
mutations <-unique(dat$mutID) # 79972
setDTthreads(1)

DF <- data.frame(expr="",N=NA,time=NA,sd=NA,se=NA,ci=NA,size=NA,stringsAsFactors = F)
DF <- DF[-1,]

summarizeMutations <- function(mut){
  # Written with data.table
  var1 <- dat[mutID == mut, .N , .(Primary.site)]
  setorder(var1, -N)
  return(list(
    mut,
    sum(var1[, .(N)]),
    stringi::stri_paste(var1$Primary.site, ":",
                        var1$N,
                        collapse = ',')
  ))
}

for(i in c(3000, 30000)){
  bmark <- microbenchmark("for" = {
    # initialize a list & counter with expected size. 
    # This speeds things considerably. For certain tasks this might not be possible.
    results_1 <- vector(mode = "list", length = i*3)
    counter <- 1 
    
    for(mut in mutations[1:i]){
      results_1[seq(counter, counter + 2)] <- summarizeMutations(mut)
      counter <- counter + 3
    }
    results_1 <- unlist(results_1, use.names = F)
  },
  "foreach" = {
    myCluster <- makeCluster(4,
                             type = "FORK",
                             useXDR = F,
                             .combine = cbind)
  print(myCluster)
  registerDoParallel(myCluster)
  results_2 <- foreach(mut = mutations[1:i], .combine = cbind, .inorder = F) %dopar% {
    return(summarizeMutations(mut))
  }
  stopCluster(myCluster)
  results_2 <- unlist(results_2, use.names = F)
  },
  "lapply"= {
    results_3 <-
      unlist(lapply(X = mutations[1:i], FUN = summarizeMutations), use.names = F)
    },
  "mclapply"={
    results_4 <- unlist(
      parallel::mclapply(
        X = mutations[1:i],
        FUN = summarizeMutations,
        mc.cores = parallel::detectCores()
      ),
      use.names = F
    )
  }, times = 5)
  # saveRDS(bmark,file = paste0("./bmark/bmark_",test.name,"_",i,".RDS"))
  results <- summarySE(bmark,measurevar = "time",groupvars = "expr",statistic = "mean")
  results$size <- rep(i,length(results[,1]))
  DF <- rbind.data.frame(DF,results)
  rm(results,bmark)
}
# saveRDS(DF,file = paste0("./results/results_",test.name,".RDS"))
rm(list=ls())
gc()

# ---------- parallel saveRDS.gz()
source("https://raw.githubusercontent.com/dchakro/shared_Rscripts/master/summarySE.R")
source('https://gist.githubusercontent.com/dchakro/8b1e97ba6853563dd0bb5b7be2317692/raw/parallelRDS.R')
DT <- readRDS.gz("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/COSMIC_v92_R_DT/allCodingMutations.RDS")

test.name <- "parallel_saveRDS"
size <- c("100", "10000", "1000000")

DF <- data.frame(expr="",N=NA,time=NA,sd=NA,se=NA,ci=NA,size=NA,stringsAsFactors = F)
DF <- DF[-1,] 

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
  saveRDS(bmark,file = paste0("./bmark/bmark_",test.name,"_",i,".RDS"))
  results <- summarySE(bmark,measurevar = "time",groupvars = "expr",statistic = "mean")
  results$size <- rep(i,length(results[,1]))
  DF <- rbind.data.frame(DF,results)
  rm(results,bmark)
}
saveRDS(DF,file = paste0("./results/results_",test.name,".RDS"))
rm(list=ls())
gc()


# ---------- parallel readRDS.gz()
source("https://raw.githubusercontent.com/dchakro/shared_Rscripts/master/summarySE.R")
source('https://gist.githubusercontent.com/dchakro/8b1e97ba6853563dd0bb5b7be2317692/raw/parallelRDS.R')
FILE <- "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/COSMIC_v92_R_DT/CountStatsRAW.RDS" # 51 MB RDS file

test.name <- "parallel_readRDS"

DF <- data.frame(expr="",N=NA,time=NA,sd=NA,se=NA,ci=NA,size=NA,stringsAsFactors = F)
DF <- DF[-1,] 

bmark <- microbenchmark(
  "base" = {
    base::readRDS(file = FILE)
  },
  "parallel" = {
    readRDS.gz(file = FILE, threads = parallel::detectCores())
  }
  ,times = 10)
gc()
saveRDS(bmark,file = paste0("./bmark/bmark_",test.name,".RDS"))
results <- summarySE(bmark,measurevar = "time",groupvars = "expr",statistic = "mean")
results$size <- rep("51 MB", length(results[, 1]))
DF <- rbind.data.frame(DF, results)
rm(results, bmark)
saveRDS(DF, file = paste0("./results/results_", test.name, ".RDS"))
rm(list = ls())
gc()

