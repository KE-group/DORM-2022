setwd("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Data/adding Targ to genomewide")
library(data.table)
#-----------------------
# Reading data
rm(list=ls()); gc()

df_gw <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/genome_wide/20220117_results/20220117.FrequencyByMutation.RDS")
df_gw <- df_gw[,1:4]
df_gw[, MutID:=paste0(Gene,"=",Mutation)]

df_full <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/targeted_and_wgs/20220606_results/20220606.FrequencyByMutation.RDS")
df_full <- df_full[,1:4]
df_full[, MutID:=paste0(Gene,"=",Mutation)]

# -----------
# Generating data for Plotting
N=5000
selection <- data.frame(MutID=df_gw$MutID[1:N])
selection$gw <- df_gw$counts[match(x = selection$MutID, table = df_gw$MutID)]
selection$full <- df_full$counts[match(x = selection$MutID, table = df_full$MutID)]
selection <- setDT(selection)

# # Removing genome-wide data from full and renaming to targeted seq
selection[,targ := full - gw]

setnames(selection,"targ","targ_count")
setnames(selection,"full","full_count")
setnames(selection,"gw","gw_count")

# Converting to % of samples
# GW 36224; FULL 328017 ; TARG 364241
selection[, gw := ((gw_count / 36224) * 100)]
selection[, targ := ((targ_count / 328017) * 100)]
selection[, full := ((full_count / 364241) * 100)]

#-----------
#  Plotting
library(ggplot2)
source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type = "L",legend = F)
ggplot(data = selection, aes(x=gw, y=targ))+
  geom_point(alpha=1, 
             shape = 21,
             stroke = 0.2,
             color = "#000000",
             fill = "#1280C3",
             aes(size=gw_count))+
  xlab("MAF % (genome-wide sequencing)")+
  ylab("MAF % (in targeted sequencing)")+
  geom_abline(slope = 1,intercept = 0)+
  scale_x_continuous(expand=c(0,0),limits = c(0,16))+
  scale_y_continuous(expand=c(0,0),limits = c(0,16))+
  scale_size(range=c(1,10))+
  ggtitle(paste0("Top ",N, " Mutations"))+
  customtheme+
  coord_cartesian(clip = "on")

ggsave(
  filename = paste0("gg_Scatter_top", N, ".pdf"),
  height = 5,
  width = 5
)

library(plotly)
p <- plot_ly(
  selection,
  alpha = 0.8,
  x =  ~ gw,
  y =  ~ targ,
  size = ~gw_count,
  sizes = c(8, 50),
  colors = "#1280C3",
  mode = "text",
  text =  ~ paste("Cell Line: ", gsub("="," ", MutID),
                  "\nCount (T) = ", targ_count,
                  "\nCount (G) = ", gw_count)) %>%
  layout(shapes = list(list(
    type = "line",
    x0 = 0, 
    x1 = ~max(selection$gw, selection$targ), 
    xref = "x",
    y0 = 0, 
    y1 = ~max(selection$gw, selection$targ),
    yref = "y",
    line = list(color = "black", dash="dot")
  ))) %>%
  add_markers(marker=list(opacity = 0.5, sizemode = 'diameter')) %>%
  layout(
    title = paste0("<b>Top ",N, " Mutations</b>"),
    xaxis = list(title = paste0("<b>MAF in genome-wide screens (%)</b>")),
    yaxis = list(title = paste0("<b>MAF in complete data (%)</b>"))
  )
# hide_colorbar(p)
htmlwidgets::saveWidget(
  widget = hide_colorbar(p),
  file = paste0("gg_Scatter_top", N, ".html"),
  selfcontained = T,
  libdir = NULL)

# selection[, MutID := gsub("=", " ", MutID)]
