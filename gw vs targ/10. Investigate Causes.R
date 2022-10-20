setwd("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Data/adding Targ to genomewide")
library(data.table)
rm(list=ls()); gc()

# ##------- Setting up data
# 
# source('https://gist.githubusercontent.com/dchakro/8b1e97ba6853563dd0bb5b7be2317692/raw/parallelRDS.R')
# # Unused --------------
# # gwDT <- readRDS.gz("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/genome_wide/1_AllsamplesMinInfo.RDS")
# #
# # gwDT <- gwDT[!grepl("_ENST", Gene.name, fixed = T), ]
# # gwDT <- gwDT[!stringi::stri_detect_regex(str = Mutation.AA, pattern = "^p.\\?"), ]
# # gwDT <- gwDT[!grepl(pattern = '=', x = Mutation.AA, fixed = T), ]
# # gwDT[, Mutation.AA := stringi::stri_replace_first_fixed(str = Mutation.AA,
# #                                                     pattern = "p.",
# #                                                     replacement = "")]
# # gwDT[, Mutation.AA := stringi::stri_replace_all_regex(str = Mutation.AA,
# #                                                   pattern = "\\*",
# #                                                   replacement = "X")]
# 
# fullDT <- readRDS.gz("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/targeted_and_wgs/1_AllsamplesMinInfo.RDS")
# 
# fullDT <- fullDT[!grepl("_ENST", Gene.name, fixed = T), ]
# fullDT <- fullDT[!stringi::stri_detect_regex(str = Mutation.AA, pattern = "^p.\\?"), ]
# fullDT <- fullDT[!grepl(pattern = '=', x = Mutation.AA, fixed = T), ]
# fullDT[, Mutation.AA := stringi::stri_replace_first_fixed(str = Mutation.AA,
#                                                         pattern = "p.",
#                                                         replacement = "")]
# fullDT[, Mutation.AA := stringi::stri_replace_all_regex(str = Mutation.AA,
#                                                       pattern = "\\*",
#                                                       replacement = "X")]
# 
# PMID_info <- unique(na.exclude(fullDT[,.(Sample.name, PMID, Gene.name)]))
# # Calculating the scale of study i.e. number of genes assayed in the individual studies
# number_of_genes <- PMID_info[,uniqueN (.SD), 
#                              by=PMID, 
#                              .SDcols = c("PMID", "Gene.name")]
# rm(PMID_info);gc()
# 
# # Calculate the number of samples in the study (with atleast 1 mutation)
# tmp <- unique(na.exclude(fullDT[,.(Sample.name,PMID,Primary.site,Histology)]))
# 
# sampleCount_in_study <- tmp[, .(count = .N) , by = PMID]
# sampleCount_by_histology <- tmp[, .(count = .N) , by = c("PMID","Primary.site","Histology")]
# rm(tmp);gc()
# 
# ##--- Writing Data
# con <- pipe("pigz -p4 > 20221020_Data.GW.Targ.gz", "wb")
# save(fullDT,
#      sampleCount_in_study,
#      sampleCount_by_histology,
#      number_of_genes,
#      file = con)
# close(con);rm(con)


## -- Reading Data
con <- pipe("pigz -dkc -p4 20221020_Data.GW.Targ.gz","rb")
load(file = con)
close(con);rm(con)

#---- Initializing essential variables
library(ggplot2)
library(patchwork)
source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type = "L",legend = T)
dir.create("investigate/",showWarnings = F)

#---------------
# library(ggplot2)
# source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
# customtheme <- DC_theme_generator(type = "L")
# 
# ggplot(number_of_genes,aes(x=1,y=V1))+
#   geom_violin()+
#   scale_y_log10()+
#   ylab("Number of genes\nassayed in study")+
#   customtheme

#------------------
##---- Defining Functions

# # Calculate the number of samples in the study (with atleast 1 mutation)
countSamplesForPMID <- function(gene, mutation, primary.site, histology = NULL,PubMedID = NULL){
  if(is.null(histology)){
    return(nrow(unique(x = fullDT[PMID == PubMedID & Gene.name==gene & primary.site==primary.site, "Sample.name"])))
  } else {
    return(nrow(unique(x = fullDT[PMID == PubMedID & Gene.name==gene & primary.site==primary.site & Histology == histology, "Sample.name"])))
  }
}

plotStudyStats <- function(gene, mutation, primary.site, histology = NULL){
  if(is.null(histology)){
    # Selecting data related to the input Query gene/mutation/tissue combo
    subFull <- fullDT[Gene.name==gene & Mutation.AA==mutation & Primary.site == primary.site,]
    
  } else {
    # Selecting data related to the input Query gene/mutation/tissue combo
    subFull <- fullDT[Gene.name==gene & Mutation.AA==mutation & Primary.site == primary.site & Histology == histology,]
  }
      
    # Calculating number of samples with the query mutation in various studies
    alterationFreq <- subFull[,.N, .(PMID, Genomewide.screen)]
    rm(subFull);gc()
    alterationFreq <- alterationFreq[!is.na(PMID),]
    idx <- match(x = alterationFreq$PMID,table = sampleCount_in_study$PubMedID)
    alterationFreq[, percent:= N / sampleCount_in_study$N[idx]]
    
    idx <- match(x = alterationFreq$PMID,table = number_of_genes$PMID)
    alterationFreq[, studyScale:= number_of_genes$V1[idx]]
    
    setorder(alterationFreq,-N)
    set.seed(2022)
    
    print(head(alterationFreq,20))
    # # Scatterplot -  Percent altered vs PMID (size of dot = study size)
    # ggplot(alterationFreq,aes(x=sample(PMID),
    #                           y=percent,
    #                           size=N,
    #                           color=Genomewide.screen))+
    #   geom_point(alpha = 0.5,shape=16)+
    #   xlab("PubMed ID")+
    #   ggtitle(paste0(gene,
    #                  " ",
    #                  mutation,
    #                  " - ",
    #                  primary.site,
    #                  ifelse(is.null(histology), yes = "", no = paste0(" / ",histology))))+
    #   ylab(paste0("Percentage of samples in dataset\nwith ",gene," ",mutation))+
    #   scale_y_continuous(expand = c(0,0),labels = scales::percent,limits = c(0,1))+
    #   scale_x_continuous(expand = c(0,0))+
    #   scale_size(range = c(0.5,10))+
    #   scale_color_manual(values = c("#1380A1","#990000"))+
    #   customtheme+
    #   theme(axis.ticks.x = element_blank(),
    #         axis.text.x = element_blank())+
    #   coord_cartesian(clip="off")
    
    ## Scatterplot -  percent vs study size
    # ggplot(alterationFreq,aes(x=N,
    #                           y=percent,
    #                           color=Genomewide.screen))+
    #   geom_point(alpha = 0.5,shape=16,size=2)+
    #   xlab("Study size (n)")+
    #   ggtitle(paste0(gene,
    #                  " ",
    #                  mutation,
    #                  " - ",
    #                  primary.site,
    #                  ifelse(is.null(histology), yes = "", no = paste0(" / ",histology))))+
    #   ylab(paste0("Percentage of samples in dataset\nwith ",gene," ",mutation))+
    #   scale_y_continuous(expand = c(0,0),labels = scales::percent,limits = c(0,1))+
    #   scale_x_continuous(expand = c(0,0))+
    #   scale_size(range = c(0.5,10))+
    #   scale_color_manual(values = c("#1380A1","#990000"))+
    #   customtheme+
    #  coord_cartesian(clip="off")
    
    studyScale_breaks <-
      c(0, 1, 10, 30, 100, 1000, max(alterationFreq$studyScale, na.rm = T))

    bar_PopulationFreq <- ggplot(
      data = alterationFreq,
      aes(
        x = reorder(PMID, -N),
        y = percent,
        fill = studyScale
      )
    )+
      geom_col(position = "dodge",
               color=NA)+
      customtheme+
      theme(axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank())+
      scale_y_continuous(labels = scales::percent,expand = c(0,0))+
      ylab(paste0("Altered samples (%)\n(",gene," ",mutation,")"))+
      binned_scale(aesthetics = "fill",
                   scale_name = "stepsn",
                   palette = function(x) viridis::turbo(length(studyScale_breaks), direction = -1),
                   breaks = studyScale_breaks,
                   guide = "colorsteps")+
      ggtitle(paste0(
        gene,
        " ",
        mutation,
        " - ",
        primary.site,
        ifelse(
          is.null(histology),
          yes = "",
          no = paste0(" / ", histology)
        )
      ))

    bar_Size <- ggplot(data = alterationFreq,
                       aes(x = reorder(PMID, -N),
                           y = N,
                           fill = studyScale))+
      geom_col(position = "dodge",
               color=NA)+
      customtheme+
      theme(axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.line.x = element_blank())+
      scale_y_reverse(expand = c(0,0))+
      ylab("Study size (n)")+
      geom_hline(yintercept = 0,size=1)+
      binned_scale(aesthetics = "fill",
                   scale_name = "stepsn",
                   palette = function(x) viridis::turbo(length(studyScale_breaks),direction = -1),
                   breaks = studyScale_breaks,
                   guide = "colorsteps")

    export <-  bar_PopulationFreq +
      bar_Size +
      plot_layout(nrow = 2,
                  ncol=1,
                  widths =  unit(c(28, 28), c('cm', 'cm')),
                  heights = unit(c(4, 4), c('cm', 'cm')),
                  guides = "collect")

    ggsave(
      filename = paste0(
        "investigate/",
        ifelse(
          test = is.null(histology),
          yes = paste(gene,
                      mutation,
                      primary.site,
                      sep = "_"),
          no = paste(gene,
                     mutation,
                     primary.site,
                     gsub("/","_",histology,fixed = T),
                     sep = "_")),
          ".pdf"),
      plot = export,
      width = 35,
      height = 10,
      units = "cm"
      )
    return(alterationFreq)
}

val <- plotStudyStats("EGFR","L858R","lung","adenocarcinoma")
plotStudyStats("EGFR","L858R","lung")
plotStudyStats("EGFR","E746_A750del","lung")
plotStudyStats("KRAS","G12C","lung")
plotStudyStats("KRAS","G12V","pancreas","ductal_carcinoma")
plotStudyStats("KRAS","G12D","large_intestine")
plotStudyStats("BRAF","V600E","skin","malignant_melanoma/NS")
plotStudyStats("PIK3CA","H1047R","breast")
plotStudyStats("PIK3CA","H1047R","breast","carcinoma/NS")
plotStudyStats("JAK2","V617F","haematopoietic_and_lymphoid_tissue")
plotStudyStats("JAK2","V617F","haematopoietic_and_lymphoid_tissue","haematopoietic_neoplasm/polycythaemia_vera")
plotStudyStats("JAK2","V617F","haematopoietic_and_lymphoid_tissue","haematopoietic_neoplasm/myelofibrosis")
plotStudyStats("IDH1","R132H","central_nervous_system")

# 
gene <- "EGFR"
mutation <- "L858R"
primary.site <- "lung"
# histology <- "adenocarcinoma"

# gene <- "BRAF"
# mutation <- "V600E"
# primary.site <- "skin"
# histology <- "malignant_melanoma/NS"

gene <- "IDH1"
mutation <- "R132H" 
primary.site <- "central_nervous_system"
histology <- NULL

# rm(list=ls()[!ls() %in% c("fullDT","sampleCount_in_study","number_of_genes")])