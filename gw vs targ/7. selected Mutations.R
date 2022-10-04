setwd("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Data/adding Targ to genomewide")
library(data.table)
rm(list=ls()); gc()
#----------------
source('https://gist.githubusercontent.com/dchakro/8b1e97ba6853563dd0bb5b7be2317692/raw/parallelRDS.R')

# Define Functions
findHistology <- function(mutant = NULL, primary.site = NULL){
  ## Description: Generate prevalence of a mutant in various histological diagnoses
  
  if(grepl(pattern = "=",x = mutant,fixed = T)){
    # selecting data for the specific mutant
    gw_mut <- gw[mutID==mutant & Primary.site== primary.site,]
    full_mut <- full[mutID==mutant & Primary.site== primary.site,]
    
    # Calculating statistics for each histology
    f_stat <- full_mut[,.N,.(Histology)]
    gw_stat <- gw_mut[,.N,.(Histology)]
    f_stat[,percentage := (N*100)/nrow(full_mut)]
    gw_stat[,percentage := (N*100)/nrow(gw_mut)]
    setorder(f_stat,-percentage)
    setorder(gw_stat,-percentage)
    
    # Calculating "Out of total" stats
    
    # clearing memory
    rm(gw_mut,full_mut);gc()
    
    # Returning value
    return(list(FULL=f_stat,GenomeWide=gw_stat))
  } else {
    message(paste("Incorrect format for the parameter mutant.\n Use format \"GENE=AAchange\" like findHistology(mutant = \"BRAF=V600E\") \n Received input:",mutant))
  }
}

gather_saveData <- function(mutant = NULL, primary.site = NULL) {
  ## Description: uses findHistology() to create summary stats for full & genome-wide data
  ## Input: mutant and tissue
  ## Output: writes the tables as files
  ##         N = number of mutations in that Primary.site with that histology
  
  if(grepl(pattern = "=",x = mutant,fixed = T)){
    data_of_mutant_for_primary_site <- findHistology(mutant = mutant, primary.site = primary.site)
    site_full <- full[Primary.site==primary.site, .N, .(Histology)]
    site_gw <- gw[Primary.site==primary.site, .N, .(Histology)]
    
    DT <- data_of_mutant_for_primary_site[["FULL"]]
    site_full[, FREQ := DT$N[match(x = Histology,table = DT$Histology)]]
    site_full[,percentage := (FREQ*100)/N]
    setorder(site_full,-FREQ,na.last = T)
    setnames(x = site_full,old = "FREQ",new = gsub("=","_",mutant))
    
    DT <- data_of_mutant_for_primary_site[["GenomeWide"]]
    site_gw[, FREQ := DT$N[match(x = Histology,table = DT$Histology)]]
    site_gw[,percentage := (FREQ*100)/N]
    setorder(site_gw,-FREQ,na.last = T)
    setnames(x = site_gw,old = "FREQ",new = gsub("=","_",mutant))
    
    write.table(
      x = site_full,
      file = paste0("mutants/",primary.site,"_",gsub("=","_",mutant),"_full.tsv"),
      sep = "\t",
      quote = F,
      row.names = F,
      col.names = T
    )
    
    write.table(
      x = site_gw,
      file = paste0("mutants/",primary.site,"_",gsub("=","_",mutant),"_gw.tsv"),
      sep = "\t",
      quote = F,
      row.names = F,
      col.names = T
    )
  } else {
    message(paste("Incorrect format for the parameter mutant.\n Use format \"GENE=AAchange\" like findHistology(mutant = \"BRAF=V600E\") \n Received input:",mutant))
  }
}


#------ Set up data
gw <- readRDS.gz("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/genome_wide/4.COSMIC.all.coding.Mutatations.RDS")
gw[,Mutation.AA := stringi::stri_replace_first_fixed(str = Mutation.AA,
                                                     pattern = "p.",
                                                     replacement = "")]
gw[,Mutation.AA := stringi::stri_replace_first_fixed(str = Mutation.AA,
                                                     pattern = "\\*",
                                                     replacement = "X")]
gw[,Gene.name:= as.character(Gene.name)]
gw[, mutID := paste(Gene.name, Mutation.AA, sep = "=")]

full <- readRDS.gz("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/targeted_and_wgs/4.COSMIC.all.coding.Mutatations.RDS")
full[,Mutation.AA := stringi::stri_replace_first_fixed(str = Mutation.AA,
                                                       pattern = "p.",
                                                       replacement = "")]
full[,Mutation.AA := stringi::stri_replace_first_fixed(str = Mutation.AA,
                                                       pattern = "\\*",
                                                       replacement = "X")]
full[,Gene.name:= as.character(Gene.name)]
full[, mutID := paste(Gene.name, Mutation.AA, sep = "=")]

#--------
# Mutations to analyze
## JAK2=V617F - haematopoietic_and_lymphoid_tissue
## EGFR=L858R - lung
## EGFR=E746_A750del - lung
## KRAS=G12V - pancreas
## BRAF=V600E - thyroid

#--------------
rm(list=ls()[!ls() %in% c("full","gw","findHistology" , "gather_saveData")])
gather_saveData(mutant = "JAK2=V617F", primary.site = "haematopoietic_and_lymphoid_tissue")
gather_saveData(mutant = "EGFR=L858R", primary.site = "lung")
gather_saveData(mutant = "EGFR=E746_A750del", primary.site = "lung")
gather_saveData(mutant = "KRAS=G12V", primary.site = "pancreas")
gather_saveData(mutant = "BRAF=V600E", primary.site = "thyroid")

