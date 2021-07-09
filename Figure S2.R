library(ggplot2)
library(data.table)
rm(list=ls());gc()
# ---> Figure S2 <----

source("https://gist.githubusercontent.com/dchakro/8b1e97ba6853563dd0bb5b7be2317692/raw/parallelRDS.R")
muts <- readRDS.gz("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v94/Full_Database/4.COSMIC.all.coding.Mutatations.RDS")
rm(loadRDS,readRDS.gz,writeRDS,saveRDS.gz)

muts[, Mutation.AA := stringi::stri_replace_first_fixed(str = Mutation.AA,
                                                        pattern = "p.",
                                                        replacement = "")]

muts[, Mutation.AA := stringi::stri_replace_all_regex(str = Mutation.AA,
                                                      pattern = "\\*",
                                                      replacement = "X")]

muts[, mutID := paste(Gene.name,Mutation.AA,sep="=")]

output <- muts[,.(.N), by = mutID]
data.table::setnames(output, "N", "freq")
singleOccurances <- output[freq==1, .(mutID)]
rm(output)

'%nin%' <- Negate('%in%')
dataDF <- muts[, .(Sample.name, mutID, Primary.site)]
dataDF <- dataDF[mutID %nin% singleOccurances$mutID, ]
rm("%nin%",singleOccurances, muts)
gc()

# ---> Figure S2 A : No. of samples with recurrent mutations by tissue type <----
rm(list = ls()[! ls() %in% c("dataDF", "DC_theme_generator")])
SampleDF <- unique(dataDF[,.(Sample.name, Primary.site)])
Stats <- SampleDF[,.(.N), by = Primary.site]
colnames(Stats) <- c("tissue","count")
setorder(Stats, -count)
Stats[, tissue := gsub("_", " ", tissue)]

source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type = 'L',
                                  legend = 'F',
                                  ticks = 'out',
                                  x.axis.angle = 45,
                                  hjust = 1,
                                  vjust = 1,
                                  fontsize.cex = 1.2,
                                  ax.fontstyle = "italic")
options(scipen=100000)
ggplot(data = Stats, aes(y=count,
                         x=reorder(tissue, -count)))+
  geom_col(fill="#000000",
           width=0.75)+
  xlab("Tissue of origin of cancer")+
  ylab("Number of samples\nwith recurrent mutations")+
  customtheme +
  scale_y_continuous(expand = c(0, 0))

ggsave(
  "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Figures/panels from R/S2 A.pdf",
  width = 8,
  height = 5,
  device = cairo_pdf
)


# ---> Figure S2 B : No. of Recurrent mutations by cancer type <----
rm(list = ls()[! ls() %in% c("dataDF", "DC_theme_generator")])

Stats <- dataDF[,.(.N), by = Primary.site]
colnames(Stats) <- c("tissue","count")
setorder(Stats, -count)
Stats[, tissue := gsub("_", " ", tissue)]

source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type = 'L',
                                  legend = 'F',
                                  ticks = 'out',
                                  x.axis.angle = 45,
                                  hjust = 1,
                                  vjust = 1,
                                  fontsize.cex = 1.2,
                                  ax.fontstyle = "italic")
options(scipen=100000)
ggplot(data = Stats, aes(y=count,
                         x=reorder(tissue, -count)))+
  geom_col(fill="#000000",
           width=0.75)+
  customtheme+
  scale_y_continuous(expand = c(0,0))+
  xlab("Tissue of origin of cancer")+
  ylab("Number of recurrent mutations")


ggsave(
  "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Figures/panels from R/S2 B.pdf",
  width = 8,
  height = 5,
  device = cairo_pdf
)

# ---> Figure S2 C : TMB plot <----
rm(list = ls()[! ls() %in% c("dataDF", "DC_theme_generator")])

Stats <- dataDF[,.(.N), by = Sample.name]
colnames(Stats) <- c("Sample.name","count")
setorder(Stats, -count)
Sample_Tissue_Map <- unique(dataDF[,.(Sample.name, Primary.site)])
# print("Number of samples by tissue")
# Sample_Tissue_Map[, .N, .(Primary.site)]
Stats[, tissue := Sample_Tissue_Map$Primary.site[match(x = Stats$Sample.name, table = Sample_Tissue_Map$Sample.name)]]
Stats[, tissue := gsub("_", " ", tissue)]

source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type = 'L',
                                  legend = 'F',
                                  ticks = 'out',
                                  x.axis.angle = 45,
                                  hjust = 1,
                                  vjust = 1,
                                  fontsize.cex = 1.2,
                                  ax.fontstyle = "italic")
options(scipen=100000)
ggplot(data = Stats, aes(y=count,
                         x=tissue))+
  geom_jitter(alpha=0.25,
              height = 0,
              size=1.5,
              width=0.3)+
  scale_y_continuous()+
  xlab("Tissue of origin of cancer")+
  ylab("Number of mutations\nper sample")+
  customtheme 

# ggsave(
#   "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/Barplots/TMB.pdf",
#   width = 8,
#   height = 5,
#   device = cairo_pdf
# )

ggsave(
  "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Figures/panels from R/S2 C.pdf",
  width = 8,
  height = 5,
  device = cairo_pdf
)
