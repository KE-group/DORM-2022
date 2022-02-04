#-------->for vs lapply vs foreach vs  mclapply <----------

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

bmark_out_tsv <- paste0("table/bmark","_",test.name,".tsv")
file.create(bmark_out_tsv, showWarnings = F)
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