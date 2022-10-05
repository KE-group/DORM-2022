to_install <- c("patchwork")  ## Set packages to install here and run
if (unname(Sys.info()["sysname"]) == "Darwin") {
  CPU_cores <- as.integer(system(command = "sysctl -n hw.physicalcpu", intern = T))
  ## For logical cores:
  # CPU_cores <- system("sysctl -n hw.ncpu")
} else if (unname(Sys.info()["sysname"]) == "Linux") {
  CPU_cores <- as.integer(system(command = "nproc", intern = T))
}
utils::setRepositories(ind = c(1, 2, 3))
missing_packages <-
  to_install[!(to_install %in% installed.packages()[, "Package"])]
if (length(missing_packages))
  install.packages(missing_packages, type = "source", INSTALL_opts = "--byte-compile", Ncpus = CPU_cores)
rm(missing_packages, to_install, CPU_cores)

setwd("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Data/adding Targ to genomewide/mutants/")
dir.create("plots",showWarnings = F)
library(data.table)
library(patchwork)
rm(list=ls()); gc()

source("https://raw.githubusercontent.com/dchakro/shared_Rscripts/master/roundUp.R")

library(ggplot2)
source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type = "L",legend = F)

files <- as.list(setdiff(list.files(), list.dirs(recursive = FALSE, full.names = FALSE)))
#----------
saveScatterplots <- function(fileName){
  tissue <-
    gsub(
      pattern = "_tsv",
      replacement = "",
      x = paste0(
        unlist(
          stringi::stri_extract_all_regex(str = gsub("del", "DEL", fileName , fixed = T), pattern = "[[:lower:]]+[_]*"),
          use.names = F
        ),
        collapse = ""
      )
    )
  
  mutID <-
    gsub(
      pattern = "^[_]+",
      replacement = "",
      paste0(unlist(
        stringi::stri_extract_all_regex(str = gsub("del", "DEL", fileName , fixed = T), pattern = "[^a-z.]+")
      ), collapse = ""),
      fixed = F
    )
  # reading file
  myDT <- data.table::fread(file = fileName ,
        sep = "\t",
        skip = 0,
        header = T,
        stringsAsFactors = F,
        showProgress = T,
        nThread = 1)

  # Keeping rows where there are at least 10 samples in either "threshold" samples in 
  # gw or full i.e., gw+targ data. AND where there is at least 1 sample analyzed in the gw data.
  myDT <- myDT[which((!is.na(myDT[, 2])) & 
                 (myDT[, 3] > 10 | myDT[, 6] > 10)), ]
  
  colnames(myDT)[3] <- "gw_mutant_N"
  colnames(myDT)[6] <- "full_mutant_N"
  
  ax_limit <- roundUp(max(myDT$gw_percentage, myDT$full_percentage, na.rm = T))
  color_palette <-
    ifelse(test = length(unique(myDT$Histology)) > 9,
           yes = "Paired",
           no = "Set1")
  
  # ggplot(data <- myDT,aes(x=gw_percentage,
  #                         y=full_percentage,
  #                         size=full_mutant_N,
  #                         alpha=gw_mutant_N,
  #                         color=reorder(Histology,-full_mutant_N)))+
  #   geom_abline(slope = 1,intercept = 0, linetype="dotted")+
  #   geom_point()+
  #   scale_x_continuous(expand = c(0,0), limits = c(0,ax_limit))+
  #   scale_y_continuous(expand = c(0,0), limits = c(0,ax_limit))+
  #   xlab("Share in genome-wide screens (%)")+
  #   ylab("Share in complete data (%)")+
  #   scale_colour_viridis_d(option = "plasma")+
  #   customtheme
  
#   ggplot(data <- myDT,aes(x=gw_percentage,
#                           y=full_percentage,
#                           color=full_mutant_N,
#                           size=gw_mutant_N,
#                           shape=reorder(Histology,-full_mutant_N)))+
#     geom_abline(slope = 1,intercept = 0, linetype="dotted")+
#     geom_point()+
#     scale_x_continuous(expand = c(0,0), limits = c(0,ax_limit))+
#     scale_y_continuous(expand = c(0,0), limits = c(0,ax_limit))+
#     xlab("Share in genome-wide screens (%)")+
#     ylab("Share in complete data (%)")+
# #    scale_colour_viridis_c(option = "turbo")+
#     scale_colour_viridis_c(option = "plasma","N (full)")+
#     scale_shape(solid = TRUE, "Histology")+
#     scale_size("N (gw)")+
#     customtheme
  
  myDT[is.na(myDT)]=0
  main_plot <- ggplot(data <- myDT,aes(x=gw_percentage,
                          y=full_percentage,
                          size=full_mutant_N,
                          alpha=gw_mutant_N,
                          color=reorder(Histology,-full_mutant_N)))+
    geom_abline(slope = 1,intercept = 0, linetype="dotted")+
    geom_point()+
    scale_x_continuous(expand = c(0,0), limits = c(0,ax_limit))+
    scale_y_continuous(expand = c(0,0), limits = c(0,ax_limit))+
    xlab("Share in genome-wide screens (%)")+
    ylab("Share in complete data (%)")+
    ggtitle(gsub("_"," ", paste0(mutID," in ", tissue)))+
    # scale_color_brewer(palette = color_palette,"Histology")+
    scale_colour_viridis_d(option = "turbo", "Histology", direction = -1)+
    scale_size("N (full) / size")+
    scale_alpha("N (gw) / opacity",range = c(0.3, 1))
  
  scatterPlot=main_plot+customtheme
  legend=ggpubr::as_ggplot(ggpubr::get_legend(main_plot))
  
  # Fix plot size and legend takes rest of the space.
  # Plot fixed size, legend dynamic size.
  export <-  scatterPlot + 
    legend + 
    plot_layout(nrow = 2, ncol=1,
                heights = unit(c(8, 1), c('cm', 'null')))
  ggsave(
    plot = export,
    filename = paste0(
      "plots/",
      gsub(
        pattern = ".tsv",
        replacement = ".pdf",
        x = fileName
      )
    ),
    device = cairo_pdf,
    width = 5,
    height = 11
  )
  
  
  # ggsave(
  #   plot = gridExtra::marrangeGrob(
  #     myplots,
  #     layout_matrix = matrix(
  #       data = 1:2,
  #       nrow = 2,
  #       ncol = 1,
  #       byrow = T
  #     ),
  #     as.table = F, top = NA
  #   ),
  #   filename = paste0(
  #     "plots/",
  #     gsub(
  #       pattern = ".tsv",
  #       replacement = ".pdf",
  #       x = fileName
  #     )
  #   ),
  #   device = cairo_pdf,
  #   width = 5,
  #   height = 10 + (0.2*(length(unique(myDT$Histology))))
  # )
  
  # print(paste(mutID,tissue,sep=" / "))
  return(NULL)
}

# parallel::mclapply(
#   X = files,
#   FUN = function(X)
#     saveScatterplots(X),
#   mc.cores = parallel::detectCores()
# )


lapply(
  files,
  FUN = function(X)
    saveScatterplots(X)
)

