library(ggplot2)
# ---> Hotspot Mutations: Barplot top 100 mutations <----

rm(list=ls())
resDF <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v92/Full_Database/FrequencyByMutation.RDS")
resDF <- resDF[order(-resDF$counts,resDF$Gene,resDF$Mutation,decreasing = F),]
resDF$mutsID <- paste(resDF$Gene,resDF$Mutation,sep="_")

theme_plot=theme(axis.line = element_line(colour = "black",size=0.5),panel.border = element_blank(),panel.background=element_blank(),panel.grid.major=element_blank(),axis.text.y= element_text(size = rel(1.4),color="black",margin=unit(c(0.3,0.3,0.3,0.3), "cm")),legend.key= element_rect(fill=NA,colour = NA), axis.ticks.y =element_line(colour = "black"), axis.ticks.x = element_blank(),axis.text.x = element_blank(), legend.position="none",text=element_text(family="serif"),axis.ticks.length =unit(0.2, "cm"),axis.title.y = element_text(size=rel(1.5),face="italic"),axis.title.x=element_text(size=rel(1.5),face="italic"))

n <- 100
topN <- resDF[1:n,]

ggplot(data=topN,aes(x=reorder(mutsID,-counts),y=counts))+geom_col(fill="#00abea",color=NA)+ylab("Number of somatic mutations")+xlab(paste0("Mutations (n=",length(topN$mutsID),")"))+scale_y_continuous(limits = c(0,max(topN$counts)),breaks = seq(0,1500,400),expand = c(0, 0))+theme_plot+geom_hline(yintercept = 60,linetype="dashed")

# ggsave("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/Barplots/HotspotMutations.svg",width = 5,height = 5)

n <- 100
topN <- resDF[1:n,1:3]
head(sort(table(topN$Gene),decreasing = T))
DF <- aggregate(counts ~ Gene,data=topN,FUN = function(X) sum(X))

gene_census <- vroom::vroom('/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC Cancer Gene census/Census All.tsv')
gene_census$`Role in Cancer`[gene_census$`Gene Symbol` == "TP53"] <- "TSG"

topN$Role <- gene_census$`Role in Cancer`[match(topN$Gene,gene_census$`Gene Symbol`)]
DF$Role <- gene_census$`Role in Cancer`[match(DF$Gene,gene_census$`Gene Symbol`)]

topN$Role <- gsub(", fusion","",topN$Role,fixed = T)
DF$Role <- gsub(", fusion","",DF$Role,fixed = T)
table(topN$Role)

aggregate(counts ~ Role,data=topN,FUN = function(X) sum(X))
write.table(x = topN[1:20,1:4],file ="~/Desktop/tmp.csv",sep=",",row.names = F,col.names = T)

# ---> Hotspot Residues: Barplot top 100 mutations <----

rm(list=ls())
resDF <- readRDS("/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC/v92/Full_Database/20201007.FrequencyByResidue.RDS")
resDF <- resDF[order(-resDF$counts,resDF$Gene,resDF$Mutation,decreasing = F),]
resDF$mutsID <- paste(resDF$Gene,resDF$Mutation,sep="_")

theme_plot=theme(axis.line = element_line(colour = "black",size=0.5),panel.border = element_blank(),panel.background=element_blank(),panel.grid.major=element_blank(),axis.text.y= element_text(size = rel(1.4),color="black",margin=unit(c(0.3,0.3,0.3,0.3), "cm")),legend.key= element_rect(fill=NA,colour = NA), axis.ticks.y =element_line(colour = "black"), axis.ticks.x = element_blank(),axis.text.x = element_blank(), legend.position="none",text=element_text(family="serif"),axis.ticks.length =unit(0.2, "cm"),axis.title.y = element_text(size=rel(1.5),face="italic"),axis.title.x=element_text(size=rel(1.5),face="italic"))

n <- 100
topN <- resDF[1:n,]

ggplot(data=topN,aes(x=reorder(mutsID,-counts),y=counts))+geom_col(fill="#ea3f00",color=NA)+ylab("Number of somatic mutations")+xlab(paste0("Mutated residues (n=",length(topN$mutsID),")"))+scale_y_continuous(limits = c(0,max(topN$counts)),breaks = seq(0,2500,500),expand = c(0, 0))+theme_plot+geom_hline(yintercept = 80,linetype="dashed")

# ggsave("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/Barplots/HotspotResidues.svg",width = 5,height = 5)

n <- 100
topN <- resDF[1:n,1:3]
head(sort(table(topN$Gene),decreasing = T))
DF <- aggregate(counts ~ Gene,data=topN,FUN = function(X) sum(X))

gene_census <- vroom::vroom('/Users/deepankar/OneDrive - O365 Turun yliopisto/ExtraWorkSync/Klaus-Lab-Data/Big Data/COSMIC Cancer Gene census/Census All.tsv')
gene_census$`Role in Cancer`[gene_census$`Gene Symbol` == "TP53"] <- "TSG"

topN$Role <- gene_census$`Role in Cancer`[match(topN$Gene,gene_census$`Gene Symbol`)]
DF$Role <- gene_census$`Role in Cancer`[match(DF$Gene,gene_census$`Gene Symbol`)]

topN$Role <- gsub(", fusion","",topN$Role,fixed = T)
DF$Role <- gsub(", fusion","",DF$Role,fixed = T)
table(topN$Role)

aggregate(counts ~ Role,data=topN,FUN = function(X) sum(X))

write.table(x = topN[1:20,1:4],file ="~/Desktop/tmp.csv",sep=",",row.names = F,col.names = T)
