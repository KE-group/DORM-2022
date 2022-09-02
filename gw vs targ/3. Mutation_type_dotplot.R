library(ggplot2)
library(data.table)
rm(list=ls());gc()

# Genome-wide

resDF <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/Full_Database/20220117.FrequencyByMutation.RDS")

resDF[,mutsID := paste(Gene, Mutation, sep="_")]

#---------------------------
rm(list=ls()[!ls()%in% "resDF"])
gene_census <- data.table::fread('/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC Cancer Gene census/Census All.tsv')
setnames(gene_census,make.names(colnames(gene_census)))

gene_census[Gene.Symbol == "TP53", Role.in.Cancer := "TSG"]
gene_census[Gene.Symbol == "BRAF", Role.in.Cancer := "oncogene"]
gene_census[Gene.Symbol == "CTNNB1", Role.in.Cancer := "oncogene"]
gene_census[Gene.Symbol == "ERBB2", Role.in.Cancer := "oncogene"]


genes <- resDF[1:1000,.N,.(Gene)]
genes2 <- resDF[1:1000, .(counts = sum(counts)), by = Gene]
genes[, counts := genes2$counts[match(x = Gene, table = genes2$Gene)]]
genes[, Role := gene_census$Role.in.Cancer[match(Gene, gene_census$Gene.Symbol)]]
genes[Role=="","Role"] <- NA
genes[is.na(Role),"Role"] <- "unclassified"
genes[Role=="fusion","Role"] <- "unclassified"
genes[, Role := gsub(", fusion", "-fusion", Role)]
genes[grep(pattern = ",",x = Role,fixed = T),"Role"] <- "ambiguous"
# genes[, Role := gsub(", ", "/", Role)]


rm(genes2)
setorder(genes,counts) 
library(ggplot2)
source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type = "L",x.axis.angle = 45)


# ggplot(genes,aes(x=N, y=counts, color=Role))+geom_point(size=3)+customtheme


rang <- c("#F3B61F", "#A31621","#e12433","#419D78","#69c19e","grey40")

# ggplot(genes,aes(x=Role, y=counts, color=Role, size=N))+geom_point(alpha=0.5)+customtheme+scale_y_continuous(trans = "log10")+scale_color_manual(values=rang)

ggplot(genes,aes(x=Role, y=counts, color=Role, size=N))+
  geom_point(alpha=0.5)+
  scale_y_continuous(breaks = c(0,500,seq(1000,6000,by=1000)),
                     labels=c("0","500",c(rbind(c(""),as.character(seq(2000,6000,by=2000))))))+
  xlab("Function (COSMIC Gene Census)")+
  ylab("Cumulative population frequency")+
  scale_color_manual(values=rang)+
  customtheme


ggsave(
  "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Figures/panels from R/genesVfunction_gw_only.pdf",
  width = 5,
  height = 6,
  device = cairo_pdf
)

genes[,.(ct=sum(counts)),.(Role)]


# --------------
# targeted+genomewide
rm(list=ls());gc()
resDF <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/targeted_and_wgs/20220606.FrequencyByMutation.RDS")

resDF[,mutsID := paste(Gene, Mutation, sep="_")]

#---------------------------
rm(list=ls()[!ls()%in% "resDF"])
gene_census <- data.table::fread('/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC Cancer Gene census/Census All.tsv')
setnames(gene_census,make.names(colnames(gene_census)))

# gene_census[Gene.Symbol == "TP53",.(Role.in.Cancer)]
gene_census[Gene.Symbol == "TP53", Role.in.Cancer := "TSG"]
gene_census[Gene.Symbol == "BRAF", Role.in.Cancer := "oncogene"]
gene_census[Gene.Symbol == "CTNNB1", Role.in.Cancer := "oncogene"]
gene_census[Gene.Symbol == "ERBB2", Role.in.Cancer := "oncogene"]
gene_census[Gene.Symbol == "JAK2", Role.in.Cancer := "oncogene"]


genes <- resDF[1:1000,.N,.(Gene)]
genes2 <- resDF[1:1000, .(counts = sum(counts)), by = Gene]
genes[, counts := genes2$counts[match(x = Gene, table = genes2$Gene)]]
genes[, Role := gene_census$Role.in.Cancer[match(Gene, gene_census$Gene.Symbol)]]
genes[Role=="","Role"] <- NA
genes[is.na(Role),"Role"] <- "unclassified"
genes[Role=="fusion","Role"] <- "unclassified"
genes[, Role := gsub(", fusion", "-fusion", Role)]
genes[grep(pattern = ",",x = Role,fixed = T),"Role"] <- "ambiguous"
# genes[, Role := gsub(", ", "/", Role)]


rm(genes2)
setorder(genes,counts) 
library(ggplot2)
source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')
customtheme <- DC_theme_generator(type = "L",x.axis.angle = 45)


# ggplot(genes,aes(x=N, y=counts, color=Role))+geom_point(size=3)+customtheme


rang <- c("#F3B61F", "#A31621","#e12433","#419D78","#69c19e","grey40")

# ggplot(genes,aes(x=Role, y=counts, color=Role, size=N))+geom_point(alpha=0.5)+customtheme+scale_y_continuous(trans = "log10")+scale_color_manual(values=rang)

ggplot(genes,aes(x=Role, y=counts, color=Role, size=N))+
  geom_point(alpha=0.5)+
  scale_y_continuous(breaks = c(0,seq(5000,60000,by=5000)),
                     labels = c("0",c(rbind(c(""),as.character(seq(10000,60000,by=10000)))))) +
  xlab("Function (COSMIC Gene Census)")+
  ylab("Cumulative population frequency")+
  scale_color_manual(values=rang)+
  customtheme


ggsave(
  "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Figures/panels from R/genesVfunction_targ_and_gw.pdf",
  width = 5,
  height = 6,
  device = cairo_pdf
)

genes[,.(ct=sum(counts)),.(Role)]
