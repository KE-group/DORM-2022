#------> base vs stringi <-------

# Listing amino acid residues

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

bmark_out_tsv <- paste0("table/bmark","_",test.name,".tsv")
file.create(bmark_out_tsv, showWarnings = F)
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