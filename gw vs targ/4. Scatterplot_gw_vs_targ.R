setwd("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Data/adding Targ to genomewide")
library(data.table)
#-----------------------
# Reading data
rm(list=ls()); gc()

df_gw <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/Full_Database/20220117.FrequencyByMutation.RDS")
df_gw <- df_gw[,1:4]
df_gw[, MutID:=paste0(Gene,"=",Mutation)]

df_full <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/targeted_and_wgs/20220606.FrequencyByMutation.RDS")
df_full <- df_full[,1:4]
df_full[, MutID:=paste0(Gene,"=",Mutation)]

# -----------
# Generating data for Plotting
N=5000
selection <- data.frame(MutID=df_gw$MutID[1:N])
selection$gw <- df_gw$counts[match(x = selection$MutID, table = df_gw$MutID)]
selection$full <- df_full$counts[match(x = selection$MutID, table = df_full$MutID)]
selection <- as.data.table(selection)

# # Removing genome-wide data from full and renaming to targeted seq
# selection$full <- selection$full-ifelse(is.na(selection$gw),0,selection$gw) # creates problems with mutations like (ZNF814 A337V)

setnames(selection,"full","targ_count")
setnames(selection,"gw","gw_count")

# Converting to % of samples
selection[, gw := ((gw_count / 36224) * 100)]
# Full count = 364241
# selection[, targ := ((targ_count / 328017) * 100)]
selection[, targ := ((targ_count / 364241) * 100)]

library(MASS)
get_density <- function(x, y, ...) {
  # Description: Get density of points in 2 dimensions.
  # Source: https://slowkow.com/notes/ggplot2-color-by-density/
  # Author: Kamil Slowikowski
  
  # @param x A numeric vector.
  # @param y A numeric vector.
  # @param n Create a square n by n grid to compute density.
  
  # @return The density within each square.
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

selection[,density := get_density(x = gw, y = targ, n = N)]

source("/Users/deepankar/OneDrive - O365 Turun yliopisto/Git/Gitlab.DC/Utilities/SignedFoldChange.R")
selection[,size:=FoldChange(gw_count,targ_count)]
summary(selection$size)
selection$size[selection$size == Inf] <- 0
selection$size[is.na(selection$size)] <- 0
summary(selection$size)

#-----------
#  Plotting
library(ggplot2)
source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type = "L",legend = F)
ggplot(data = selection, aes(x=gw, y=targ, color=density))+
  geom_point(alpha=0.7,aes(size=size))+
  xlab("Share in genome-wide screens (%)")+
  ylab("Share in complete data (%)")+
  scale_color_viridis_c(option = "plasma")+
  geom_abline(slope = 1,intercept = 0)+
  scale_x_continuous(expand=c(0,0),limits = c(0,16))+
  scale_y_continuous(expand=c(0,0),limits = c(0,16))+
  scale_size(range=c(1,10))+
  ggtitle(paste0("Top ",N, " Mutations"))+
  customtheme

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
  size = ~size,
  sizes = c(8, 50),
  color = ~density,
  colors = "plasma",
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
    xaxis = list(title = paste0("<b>Share in genome-wide screens (%)</b>")),
    yaxis = list(title = paste0("<b>Share in complete data (%)</b>"))
  )
# hide_colorbar(p)
htmlwidgets::saveWidget(
  widget = hide_colorbar(p),
  file = paste0("gg_Scatter_top", N, ".html"),
  selfcontained = T,
  libdir = NULL)

# selection[, MutID := gsub("=", " ", MutID)]
