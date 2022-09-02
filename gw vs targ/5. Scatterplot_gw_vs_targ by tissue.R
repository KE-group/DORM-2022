rm(list = ls())
gc()
library(data.table)

setwd("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Data/adding Targ to genomewide")
library(data.table)
#-----------------------
# Reading data
rm(list=ls()); gc()

df_gw <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/Full_Database/20220117.FrequencyByMutation.RDS")

# nCT= sample Count By Cancer Type
nCT <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Data/COSMIC_v95_R_DT/sampleCountByCancerType.RDS")
nCT[,tissue := gsub(" ", "_", tissue)]
df_gw[, c("counts","Frequency") := NULL]

df_gw <- as.data.frame(df_gw)
for(j in 3:ncol(df_gw)){
  df_gw[,j] <- (df_gw[,j]*100)/nCT$count[nCT$tissue==colnames(df_gw)[j]]
}
df_gw <- as.data.table(df_gw)
setorder(df_gw,-skin)
df_gw[, MutID:=paste0(Gene,"=",Mutation)]

df_full <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/targeted_and_wgs/20220606.FrequencyByMutation.RDS")
# nCT= sample Count By Cancer Type
nCT <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Data/COSMIC_v95_targ_and_gw/sampleCountByCancerType.RDS")
nCT[,tissue := gsub(" ", "_", tissue)]
df_full[, c("counts","Frequency") := NULL]

df_full <- as.data.frame(df_full)
for(j in 3:ncol(df_full)){
  df_full[,j] <- (df_full[,j]*100)/nCT$count[nCT$tissue==colnames(df_full)[j]]
}
df_full <- as.data.table(df_full)
setorder(df_full,-skin)

df_full[, MutID:=paste0(Gene,"=",Mutation)]

# -----------
# Generating data for Plotting
N=1000
selection <- data.frame(MutID = df_gw$MutID[1:N])
long_gw <- melt(df_gw[1:N, ], id.vars = c("MutID", "Gene", "Mutation"))
long_full <-
  melt(df_full[match(x = selection$MutID, table = df_full$MutID)], id.vars = c("MutID", "Gene", "Mutation"))

long_full[,MutID:=paste0(Gene,"=",Mutation,"_",variable)]
long_gw[,MutID:=paste0(Gene,"=",Mutation,"_",variable)]

selection <- data.frame(MutID = unique(c(long_full$MutID,long_gw$MutID)))
selection$gw <- long_gw$value[match(x = selection$MutID, table = long_gw$MutID)]
selection$full <- long_full$value[match(x = selection$MutID, table = long_full$MutID)]
selection[is.na(selection)] <- 0

# keeping muts that are detected in either of the full vs gw setting in a particular tissue
selection <- as.data.table(selection[-which(rowSums(selection[,c(2,3)])==0),])

# # Removing genome-wide data from full and renaming to targeted seq
# selection$full <- selection$full-ifelse(is.na(selection$gw),0,selection$gw) # creates problems with mutations like (ZNF814 A337V)
setnames(selection,"full","targ")

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

#-----------
#  Plotting
library(ggplot2)
source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type = "L",legend = F)
ggplot(data = selection, aes(x=gw, y=targ, color=density))+
  geom_point(alpha=0.75,size=1)+
  xlab("Share in genome-wide screens (%)")+
  ylab("Share in targeted screens (%)")+
  scale_color_viridis_c(option = "plasma")+
  geom_abline(slope = 1,intercept = 0)+
  scale_x_continuous(expand=c(0,0),limits = c(0,100))+
  scale_y_continuous(expand=c(0,0),limits = c(0,100))+
  ggtitle(paste0("Top ",N, " Mutations"))+
  customtheme

ggsave(
  filename = paste0("gg_Scatter_tissue_top", N, ".pdf"),
  height = 5,
  width = 5
)

library(plotly)
p <- plot_ly(
  selection,
  alpha = 0.8,
  x =  ~ gw,
  y =  ~ targ,
  color = ~density,
  colors = "plasma",
  mode = "text",
  text =  ~ paste("Mutant: ", gsub("[=_]"," ", MutID),
                  "\nCount (G) = ", gw,
                  "\nCount (T) = ", targ)) %>% # write counts
  layout(shapes = list(list(
    type = "line", 
    x0 = 0, 
    x1 = ~max(selection$gw, selection$targ), 
    xref = "x",
    y0 = 0, 
    y1 = ~max(selection$gw, selection$targ),
    yref = "y",
    line = list(color = "black")
  ))) %>%
  add_markers(marker = list(size = 10)) %>%
  layout(
    title = paste0("<b>Top ",N, " Mutations</b>"),
    xaxis = list(title = paste0("<b>Share in genome-wide screens (%)</b>")),
    yaxis = list(title = paste0("<b>Share in targeted screens (%)</b>"))
  )
hide_colorbar(p)
htmlwidgets::saveWidget(
  widget = hide_colorbar(p),
  file = paste0("gg_Scatter_tissue_top", N, ".html"),
  selfcontained = T,
  libdir = NULL)

# selection[, MutID := gsub("=", " ", MutID)]
