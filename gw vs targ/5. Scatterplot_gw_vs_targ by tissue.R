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

df_gw_bak <- df_gw
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

df_full_bak <- df_full
df_full <- as.data.frame(df_full)
for(j in 3:ncol(df_full)){
  df_full[,j] <- (df_full[,j]*100)/nCT$count[nCT$tissue==colnames(df_full)[j]]
}
df_full <- as.data.table(df_full)
setorder(df_full,-skin)

df_full[, MutID:=paste0(Gene,"=",Mutation)]
setnames(nCT,"count","full")
nCT$gw <- nCT_gw$count[match(nCT$tissue,nCT_gw$tissue)]
rm(nCT_gw)

# -----------
# Generating data for Plotting
N=1000
selection <- data.frame(MutID = df_gw$MutID[1:N]) # check what this line does ????
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


# selection[,density := get_density(x = gw, y = full, n = N)]
tmp_data <- as.data.frame(stringi::stri_split_fixed(str = selection$MutID,pattern = "_",simplify = T))
tmp_data <- within(tmp_data,  tissue <- paste(V2,V3,V4,V5,V6, sep=" "))
selection$Mutant <- gsub("[=]"," ", tmp_data[,1])
selection$tissue <- gsub("\\s+"," ",tmp_data$tissue)
selection$tissue <- gsub("\\s$","",selection$tissue)
rm(tmp_data)

# levels(as.factor(selection$tissue))

#-----------
#  Plotting
# colors <- viridis::turbo(length(unique(selection$tissue)))
colors_data <- data.table::fread(file = "color_vector.tsv")
colors <- colors_data$color
names(colors) <- colors_data$tissue
  
library(ggplot2)
source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type = "L",legend = F)
ggplot(data = selection, aes(x=gw, y=full, color=tissue))+
  geom_smooth(method = "lm",formula = y ~ x, aes(group=1),se = F,na.rm = T,color="#BE0000")+
  geom_abline(slope = 1,intercept = 0, linetype="dotted")+
  geom_point(alpha=0.75,size=1.5)+
  xlab("Share in genome-wide screens (%)")+
  ylab("Share in complete data (%)")+
  scale_color_manual(values=colors)+
  scale_x_continuous(expand=c(0,0),limits = c(0,100))+
  scale_y_continuous(expand=c(0,0),limits = c(0,100))+
  ggtitle(paste0("Top ",N, " Mutations"))+
  customtheme

ggsave(
  filename = paste0("gg_Scatter_tissue_top", N, ".pdf"),
  height = 5,
  width = 5
)

# #---- Plotly figure
df_gw_bak[, MutID:=paste0(Gene,"=",Mutation)]
df_full_bak[, MutID:=paste0(Gene,"=",Mutation)]

findCount_gw <- function(mutationID = NULL,
                         sampleSet = NULL) {
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

# Single-core
# unlist(lapply(
#   selection$MutID,
#   FUN = function(X)
#     findCount_gw(X, "full-data")
# ))

# Parallel processing
selection$gw_count <- unlist(
  parallel::mclapply(
    X = selection$MutID,
    FUN = function(X)
      findCount_gw(X, "genome-wide"),
    mc.cores = parallel::detectCores()
  )
)

selection$full_count <-
  unlist(
    parallel::mclapply(
      X = selection$MutID,
      FUN = function(X)
        findCount_gw(X, "full-data"),
      mc.cores = parallel::detectCores()
    )
  )

# paste0(formatC(signif(selection$gw,digits=3), digits=3,format="fg")," %")

library(plotly)
p <- plot_ly(
  data=selection,
  alpha = 0.8,
  x =  ~ gw,
  y =  ~ full,
  color = ~tissue,
  colors = colors,
  mode = "text",
  text =  ~ paste("Mutant: ", paste0(Mutant," / ",tissue),
                  "\n(Full): ",
                  paste0(formatC(signif(full, digits = 3), digits = 3, format = "fg"),
                         " % (n=",
                         full_count,
                         ")"),
                  "\n (GW): ", 
                  paste0(formatC(signif(gw, digits = 3), digits = 3, format = "fg"), 
                         " % (n=", 
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
  add_markers(marker = list(size = 10)) %>%
  layout(
    title = paste0("<b>Top ",N, " Mutations</b>"),
    xaxis = list(title = paste0("<b>Share in genome-wide screens (%)</b>")),
    yaxis = list(title = paste0("<b>Share in combined screen data (%)</b>"))
  )

# p %>% add_trace(x = ~full, y = fitted(lm(selection$gw~selection$full), mode = "lines"))
  
fv <- fitted(lm(selection$gw~selection$full))

p %>%  add_lines(
  x = selection$gw,
  y = fv,
  group = 1,
  name = "",
  line = list(
    color = "#BE0000",
    width = 0.5,
    dash = "dot"
  )
)

hide_legend(p)
htmlwidgets::saveWidget(
  widget = hide_legend(p),
  file = paste0("gg_Scatter_tissue_top", N, ".html"),
  selfcontained = T,
  libdir = NULL)

# selection[, MutID := gsub("=", " ", MutID)]
