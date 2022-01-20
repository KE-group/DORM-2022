library(ggplot2)
library(data.table)
rm(list=ls());gc()
# ---> Figure 1A : Barplot top 100 recurrent mutations <----

resDF <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/Full_Database/20220117.FrequencyByMutation.RDS")

resDF[,mutsID := paste(Gene, Mutation, sep="_")]

theme_plot=theme(axis.line = element_line(colour = "black",
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

n <- 100
topN <- resDF[1:n,]

ggplot(data=topN,aes(x=reorder(mutsID,-counts),
                     y=counts))+
  geom_col(fill="#ea3f00",
           color=NA)+
  ylab("Number of somatic mutations")+
  xlab(paste0("Top ", length(topN$mutsID), " recurrent mutations"))+
  scale_y_continuous(limits = c(0, 1600),
                     breaks = seq(0,1600,400),
                     expand = c(0, 0))+
  theme_plot+
  geom_hline(yintercept = 70,
             linetype="dashed")


ggsave(
  "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Data/Barplots/HotspotMutations.pdf",
  width = 5,
  height = 5,
  device = cairo_pdf
)
ggsave(
  "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Figures/panels from R/1A.pdf",
  width = 5,
  height = 5,
  device = cairo_pdf
)

DF <- topN[, .(counts = sum(counts)), by = Gene]
setorder(DF, -counts)

gene_census <- data.table::fread('/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC Cancer Gene census/Census All.tsv')
setnames(gene_census, make.names(colnames(gene_census)))

# gene_census[Gene.Symbol == "TP53",.(Role.in.Cancer)]
gene_census[Gene.Symbol == "TP53", Role.in.Cancer := "TSG"]

topN[, Role := gene_census[match(topN$Gene, gene_census$Gene.Symbol), .(Role.in.Cancer)]]
DF[, Role := gene_census$Role.in.Cancer[match(DF$Gene, gene_census$Gene.Symbol)]]

topN[,Role := gsub(", fusion", "", Role, fixed = T)]
DF[,Role := gsub(", fusion", "", Role, fixed = T)]
table(DF$Role)

topN[,.(sum = sum(counts)), by = Role]

write.table(x = topN[1:20,1:3],
            file ="/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Figures/panels from R/1A.table.tsv",
            sep=",",
            row.names = F,
            col.names = T)

# ---> Figure 1B Barplot top 100 recurrentlt mutated residues <----

rm(list=ls());gc()
resDF <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v95/Full_Database/20220119.FrequencyByResidue.RDS")
resDF[,mutsID := paste(Gene, Residue, sep="_")]

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

n <- 100
topN <- resDF[1:n,]

ggplot(data=topN,
       aes(x=reorder(mutsID,-counts),
           y=counts))+
  geom_col(fill="#5529C5",color=NA)+
  ylab("Number of somatic mutations")+
  xlab(paste0("Top ",length(topN$mutsID)," recurrently mutated residues"))+
  scale_y_continuous(limits = c(0,2500),
                     breaks = seq(0,2500,500),
                     expand = c(0, 0))+
  theme_plot+
  geom_hline(yintercept = 80,
             linetype="dashed")

ggsave(
  "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Data/Barplots/HotspotResidues.pdf",
  width = 5,
  height = 5,
  device = cairo_pdf
)

ggsave(
  "/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Figures/panels from R/1B.pdf",
  width = 5,
  height = 5,
  device = cairo_pdf
)

DF <- topN[, .(counts = sum(counts)), by = Gene]
setorder(DF, -counts)

gene_census <- data.table::fread('/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC Cancer Gene census/Census All.tsv')
setnames(gene_census,make.names(colnames(gene_census)))

# gene_census[Gene.Symbol == "TP53",.(Role.in.Cancer)]
gene_census[Gene.Symbol == "TP53", Role.in.Cancer := "TSG"]

topN[, Role := gene_census[match(topN$Gene, gene_census$Gene.Symbol), .(Role.in.Cancer)]]
DF[, Role := gene_census$Role.in.Cancer[match(DF$Gene, gene_census$Gene.Symbol)]]

topN[,Role := gsub(", fusion", "", Role, fixed = T)]
DF[,Role := gsub(", fusion", "", Role, fixed = T)]
table(DF$Role)

topN[,.(sum = sum(counts)), by = Role]

write.table(x = topN[1:20,1:3],
            file ="/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/DORM database/Figures/panels from R/1B.table.tsv",
            sep=",",
            row.names = F,
            col.names = T)

