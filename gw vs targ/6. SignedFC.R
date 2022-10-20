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
N=100
selection <- data.frame(MutID=df_gw$MutID[1:N])
selection$gw <- df_gw$counts[match(x = selection$MutID, table = df_gw$MutID)]
selection$full <- df_full$counts[match(x = selection$MutID, table = df_full$MutID)]
selection <- setDT(selection)

# # Removing genome-wide data from full and renaming to targeted seq
# selection$full <- selection$full-ifelse(is.na(selection$gw),0,selection$gw) # creates problems with mutations like (ZNF814 A337V)

setnames(selection,"full","targ_count")
setnames(selection,"gw","gw_count")

# Converting to % of samples
selection[, gw := ((gw_count / 36224) * 100)]
# Full count = 364241
# selection[, targ := ((targ_count / 328017) * 100)]
selection[, targ := ((targ_count / 364241) * 100)]

source("/Users/deepankar/OneDrive - O365 Turun yliopisto/Git/Gitlab.DC/Utilities/SignedFoldChange.R")
selection[,FC:=FoldChange(gw_count,targ_count)]
summary(selection$FC)
selection$FC[selection$FC == Inf] <- 0
selection$FC[is.na(selection$FC)] <- 0
summary(selection$FC)

selection$MutID <- gsub("="," ",selection$MutID)
#-----------
#  Plotting
library(ggplot2)
source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type = "L",x.axis.angle = 45,legend = F)

theme_plot <- theme(axis.line = element_line(colour = "black",
                                             size=0.5),
                    panel.border = element_blank(),
                    panel.background=element_blank(),
                    panel.grid.major=element_blank(),
                    axis.text.y= element_text(size = rel(1.4),
                                              color="black",
                                              margin=unit(c(0.3,0.3,0.3,0.3), "cm")),
                    legend.key= element_rect(fill=NA,
                                             colour = NA),
                    axis.ticks.y =element_line(colour = "black"),
                    axis.ticks.x = element_blank(),
                    axis.text.x = element_blank(),
                    legend.position="none",
                    text=element_text(family="serif"),
                    axis.ticks.length =unit(0.2, "cm"),
                    axis.title.y = element_text(size=rel(1.5),
                                                face="italic"),
                    axis.title.x=element_text(size=rel(1.5),
                                              face="italic"))
require(ggpubr)
ggdotchart(selection, x = "MutID", y = "FC",
           color = "#990000",                                # Color by groups
           sorting = "none",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           rotate = TRUE,                                # Rotate vertically
           dot.size = 2,                                 # Large dot size
           ylab = "Fold change\n(Genome-wide vs Full)",
           xlab = "Mutation",
           ggtheme = customtheme+ theme(axis.text.y= element_text(size = rel(0.8)))
)

ggsave(
  filename = paste0("gg_FoldChange_", N, ".pdf"),
  height = 10,
  width = 5
)

longDF <- reshape2::melt(selection, id.vars=c("MutID"),measure.vars=c("gw_count","targ_count"))

ggplot(data = longDF, aes(x = value, y = factor(MutID, levels = selection$MutID))) +
  geom_line(aes(group = MutID),color="grey50", linetype = "dashed") +
  geom_point(aes(color = variable))+
  scale_color_manual(values = c(gw_count="#1380A1",targ_count="#990000"))+
  customtheme+
  ylab("Mutation")+
  xlab("Number of samples")+
  theme(axis.text.y= element_text(size = rel(0.8)),
        axis.text.x= element_text(size = rel(1.2)))

ggsave(
  filename = paste0("gg_absolute_change_", N, ".pdf"),
  height = 10,
  width = 5
)

longDF <- reshape2::melt(selection, id.vars=c("MutID"),measure.vars=c("gw","targ"))

customtheme <- DC_theme_generator(type = "L",legend = F)
ggplot(data = longDF, aes(x = value, y = factor(MutID, levels = selection$MutID))) +
  geom_line(aes(group = MutID),color="grey50", linetype = "dashed") +
  geom_point(aes(color = variable))+
  scale_x_continuous(breaks=seq(0,max(longDF$value),2))+
  scale_color_manual(values = c(gw="#1380A1",targ="#990000"))+
  ylab("Mutation")+
  xlab("Share in dataset (%)")+
  customtheme+
  theme(axis.text.y= element_text(size = rel(0.8)))
    
ggsave(
  filename = paste0("gg_percentage_change_", N, ".pdf"),
  height = 10,
  width = 5
)
