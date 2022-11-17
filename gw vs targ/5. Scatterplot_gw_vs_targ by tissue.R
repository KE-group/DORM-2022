setwd("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Data/adding Targ to genomewide")
library(data.table)
#-----------------------
# Defining funcitons
rm(list=ls()); gc()

findCount_gw <- function(mutationID = NULL,
                         sampleSet = NULL) {
  ## Description: uses findHistology() to create summary stats for full & genome-wide data
  ## Input: mutation and the sampleSet
  ## Output: writes the plots to files
  
  
  tmp_data <- unlist(stringi::stri_split_fixed(str = mutationID, pattern = "_"))
  Mutant <- tmp_data[1]
  Tissue <- paste0(tmp_data[-1], collapse = "_")
  if (sampleSet == "genome-wide") {
    tryCatch(
      expr = {
        return(df_gw_bak[MutID == Mutant, get(Tissue)])
        
      },
      error = function (e) {
        message(paste("No entry found for", mutationID))
        # Actions to take
        return(NA)
      },
      finally = {
        # print("Great Success!!")
      }
    )
  } else if (sampleSet == "full-data") {
    tryCatch(
      expr = {
        return(df_full_bak[MutID == Mutant, get(Tissue)])
        
      },
      error = function (e) {
        message(paste("No entry found for", mutationID))
        # Actions to take
        return(NA)
      },
      finally = {
        # print("Great Success!!")
      }
    )
  } else {
    error("Acceptable data sets (param #2): genome-wide / full-data")
  }
}

# Reading data
# nCT = sample Count By Cancer Type
nCT_gw <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Data/COSMIC_v95_R_DT/sampleCountByCancerType.RDS")
nCT_gw[,tissue := gsub(" ", "_", tissue)]

nCT <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Data/COSMIC_v95_targ_and_gw/sampleCountByCancerType.RDS")
nCT[,tissue := gsub(" ", "_", tissue)]
setnames(nCT,"count","full")
nCT$gw <- nCT_gw$count[match(nCT$tissue,nCT_gw$tissue)]
rm(nCT_gw)

nCT[, targ:=full - ifelse(test = is.na(gw), yes = 0, no = gw)]

# Reading genome-wide data
df_gw <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/genome_wide/20220117_results/20220117.FrequencyByMutation.RDS")
df_gw[,Mutation := gsub("_", "-", Mutation)]
df_gw[, c("counts","Frequency") := NULL]

N=5000 # set number of mutations to work with

long_df_gw <- data.table::melt(data = df_gw[1:N], id.vars = c("Gene", "Mutation"))
rm(df_gw);gc()
long_df_gw <- long_df_gw[value > 10,]
long_df_gw[,MutID := paste0(Gene,"=",Mutation)]
dim(long_df_gw)
setnames(long_df_gw,"value","gw_count")
setorder(long_df_gw,-gw_count)

# Reading complete data
df_full <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/targeted_and_wgs/20220606_results/20220606.FrequencyByMutation.RDS")
df_full[,Mutation := gsub("_", "-", Mutation)]
df_full[, c("counts","Frequency") := NULL]
df_full[,MutID := paste0(Gene,"=",Mutation)]

long_df_full <-
  data.table::melt(data = df_full[match(x = unique(long_df_gw$MutID), table = df_full$MutID)], id.vars = c("Gene", "Mutation","MutID"))

# removing MutID columns as they are not needed
long_df_full[, c("MutID") := NULL]
long_df_gw[, c("MutID") := NULL]
rm(df_full);gc()

long_df_full <- long_df_full[value > 10,]
dim(long_df_full)
setnames(long_df_full,"value","full_count")
setorder(long_df_full,-full_count)


# Subtract GW counts from FULL data
long_df_gw[, MutID:=paste0(Gene,"=",Mutation,"#",as.character(variable))]
long_df_full[, MutID:=paste0(Gene,"=",Mutation,"#",as.character(variable))]

# ---------
# Generating data for Plotting
idx <- match(x = long_df_gw$MutID, table = long_df_full$MutID)
selection <- data.frame(MutID = unique(c(long_df_full$MutID,long_df_gw$MutID)))
selection$gw <- long_df_gw$gw_count[match(x = selection$MutID, table = long_df_gw$MutID)]
selection$full <- long_df_full$full_count[match(x = selection$MutID, table = long_df_full$MutID)]
rm(long_df_gw,long_df_full); gc()

# isolating Primary Site i.e., tissue from MutID
tmp_data <- as.data.frame(stringi::stri_split_fixed(str = selection$MutID,pattern = "#",simplify = T))
selection$Tissue <- tmp_data[, 2]
selection$MutID <- tmp_data[, 1]

# Separating Gene& Mutation from MutID into columns
tmp_data <- as.data.frame(stringi::stri_split_fixed(str = selection$MutID,pattern = "=",simplify = T))
selection$Gene <- tmp_data[, 1]
selection$Mutation <- tmp_data[, 2]
rm(tmp_data);gc()
setDT(selection)

selection <- selection[,c("MutID","Gene","Mutation","Tissue","gw","full")]

selection[, targ:=full - ifelse(test = is.na(gw), yes = 0, no = gw)]
setnames(x = selection,
         old = c("gw","full","targ"), 
                 new = c("gw_count","full_count","targ_count"))
# Calculating percentage of altered cases for given mutation
# in a given primary site

selection$gw <- NA
for(TISSUE in unique(selection$Tissue)){
  # Dividing the number of samples with mutation by the 
  # number of samples analyzed for that tissue type.
  selection$gw[selection$Tissue == TISSUE] <-
    selection$gw_count[selection$Tissue == TISSUE] / nCT$gw[nCT$tissue == TISSUE]
}
# multiply by 100 to convert to percentage
selection[,gw := gw*100]

selection$full <- NA
for(TISSUE in unique(selection$Tissue)){
  # Dividing the number of samples with mutation by the 
  # number of samples analyzed for that tissue type.
  selection$full[selection$Tissue == TISSUE] <-
    selection$full_count[selection$Tissue == TISSUE] / nCT$full[nCT$tissue == TISSUE]
}
# multiply by 100 to convert to percentage
selection[,full := full*100]

selection$targ <- NA
for(TISSUE in unique(selection$Tissue)){
  # Dividing the number of samples with mutation by the 
  # number of samples analyzed for that tissue type.
  selection$targ[selection$Tissue == TISSUE] <-
    selection$targ_count[selection$Tissue == TISSUE] / nCT$targ[nCT$tissue == TISSUE]
}
# multiply by 100 to convert to percentage
selection[,targ := targ*100]

# setting NAs to zero for convenience in calculation
selection[is.na(selection)] <- 0

# Retaining observations with enough n i.e. > 10
selection <- selection[full_count > 10 & gw_count > 10, ]
selection[, Tissue := gsub("_", " ", Tissue, fixed = T)]
setorder(selection, -gw_count)

write.table(
  x = selection,
  file = paste0("MAF comparison-top", N, ".csv"),
  quote = F,
  sep = ",",
  row.names = F,
  col.names = T
)

   
# -----------
# Reading color scheme
colors_data <- data.table::fread(file = "color_vector.tsv")
# colors_data[,tissue:=gsub(" ","_",tissue)]
colors <- colors_data$color
names(colors) <- colors_data$tissue

#-----------
#  Plotting
library(ggplot2)
source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type = "L",legend = F)

ggplot(data = selection, aes(x=gw, y=full, fill=Tissue))+
  geom_abline(slope = 1,intercept = 0, linetype="dotted")+
  geom_point(alpha = 1, 
             color = "#000000", 
             stroke = 0.3,
             shape = 21,
             aes(size = gw_count))+
  xlab("MAF % (genome-wide screens)")+
  ylab("MAF % (in complete data)")+
  scale_fill_manual(values = colors)+
  scale_x_continuous(expand = c(0,0), 
                     limits = c(0,100))+
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0,100))+
  scale_size(range=c(1,10))+
  coord_cartesian(clip = "off")+
  ggtitle(paste0("Most recurrent mutations (N: 1 - ",N, ")"))+
  customtheme

ggsave(
  filename = paste0("scatter_GvsF_tissue_top", N, ".pdf"),
  height = 5,
  width = 5
)

ggplot(data = selection, aes(x=gw, y=targ, fill=Tissue))+
  geom_abline(slope = 1,intercept = 0, linetype="dotted")+
  geom_point(alpha = 1, 
             color = "#000000", 
             stroke = 0.3,
             shape = 21,
             aes(size = gw_count))+
  xlab("MAF % (genome-wide sequencing)")+
  ylab("MAF % (in targeted sequencing)")+
  scale_fill_manual(values = colors)+
  scale_x_continuous(expand = c(0,0), 
                     limits = c(0,100))+
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0,100))+
  scale_size(range=c(1,10))+
  coord_cartesian(clip = "off")+
  ggtitle(paste0("Most recurrent mutations (N: 1 - ",N, ")"))+
  customtheme

ggsave(
  filename = paste0("scatter_GvsT_tissue_top", N, ".pdf"),
  height = 5,
  width = 5
)


ggplot(data = selection, aes(x=gw, y=full, fill=Tissue))+
  geom_abline(slope = 1,intercept = 0, linetype="dotted")+
  geom_point(alpha = 1, 
             color = "#000000", 
             stroke = 0.3,
             shape = 21,
             aes(size = gw_count))+
  xlab("MAF in genome-wide screens (%)")+
  ylab("MAF in complete data (%)")+
  scale_fill_manual(values = colors)+
  scale_x_continuous(expand = c(0,0), 
                     limits = c(0,25))+
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0,25))+
  scale_size(range=c(1,10))+
  coord_cartesian(clip = "off")+
  ggtitle(paste0("Most recurrent mutations (N: 1 - ",N, ")"))+
  customtheme

ggsave(
  filename = paste0("scatter_GvsF_tissue_zoom_top", N, ".pdf"),
  height = 5,
  width = 5
)
 

ggplot(data = selection, aes(x=gw, y=targ, fill=Tissue))+
  geom_abline(slope = 1,intercept = 0, linetype="dotted")+
  geom_point(alpha = 1, 
             color = "#000000", 
             stroke = 0.3,
             shape = 21,
             aes(size = gw_count))+
  xlab("MAF % (genome-wide sequencing)")+
  ylab("MAF % (in targeted sequencing)")+
  scale_fill_manual(values = colors)+
  scale_x_continuous(expand = c(0,0), 
                     limits = c(0,25))+
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0,25))+
  scale_size(range=c(1,10))+
  coord_cartesian(clip = "off")+
  ggtitle(paste0("Most recurrent mutations (N: 1 - ",N, ")"))+
  customtheme

ggsave(
  filename = paste0("scatter_GvsT_tissue_zoom_top", N, ".pdf"),
  height = 5,
  width = 5
)

# #---- Plotly figure
library(plotly)

p <- plot_ly(
  data=selection,
  size = ~gw_count,
  sizes = c(8, 40),
  x =  ~ gw,
  y =  ~ full,
  color = ~Tissue,
  colors = colors,
  mode = "text",
  text =  ~ paste("Mutant: ", paste0(Mutation," / ",Tissue),
                  "\n(Full): ",
                  paste0(formatC(signif(full, digits = 3), digits = 3, format = "fg"),
                         " % (n = ",
                         full_count,
                         ")"),
                  "\n (GW): ", 
                  paste0(formatC(signif(gw, digits = 3), digits = 3, format = "fg"), 
                         " % (n = ", 
                         gw_count, 
                         ")"))) %>% # write counts
  layout(shapes = list(list(
    type = "line", 
    x0 = 0, 
    x1 = ~max(selection$gw, selection$full), 
    xref = "x",
    y0 = 0, 
    y1 = ~max(selection$gw, selection$full),
    yref = "y",
    line = list(color = "black", dash="dot")
  ))) %>%
  add_markers(marker=list(opacity = 0.5, sizemode = 'diameter')) %>%
  layout(
    title = paste0("<b>Most recurrent mutations (N: 1 - ",N, ")</b>"),
    xaxis = list(title = paste0("<b>MAF in genome-wide screens (%)</b>")),
    yaxis = list(title = paste0("<b>MAF in combined screen data (%)</b>"))
  )

htmlwidgets::saveWidget(
  widget = hide_legend(p),
  file = paste0("scatter_GvsF_tissue_top", N, ".html"),
  selfcontained = T,
  libdir = NULL)


p <- plot_ly(
  data=selection,
  size = ~gw_count,
  sizes = c(8, 40),
  x =  ~ gw,
  y =  ~ targ,
  color = ~Tissue,
  colors = colors,
  mode = "text",
  text =  ~ paste("Mutant: ", paste0(Mutation," / ",Tissue),
                  "\n(Targ): ",
                  paste0(formatC(signif(targ, digits = 3), digits = 3, format = "fg"),
                         " % (n = ",
                         targ_count,
                         ")"),
                  "\n (GW): ", 
                  paste0(formatC(signif(gw, digits = 3), digits = 3, format = "fg"), 
                         " % (n = ", 
                         gw_count, 
                         ")"))) %>% # write counts
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
    title = paste0("<b>Most recurrent mutations (N: 1 - ",N, ")</b>"),
    xaxis = list(title = paste0("<b>MAF in genome-wide screens (%)</b>")),
    yaxis = list(title = paste0("<b>MAF in combined screen data (%)</b>"))
  )

htmlwidgets::saveWidget(
  widget = hide_legend(p),
  file = paste0("scatter_GvsT_tissue_top", N, ".html"),
  selfcontained = T,
  libdir = NULL)
