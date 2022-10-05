rm(list=ls())
setwd("/Users/deepankar/OneDrive - O365 Turun yliopisto/Klaus lab/Manuscripts/Hotspot Explorer/Data/Code optimization")
library(ggplot2)
source('https://raw.githubusercontent.com/dchakro/ggplot_themes/master/DC_theme_generator.R')

options(stringsAsFactors = F,scipen = 100000)

# ------> base_V_stringi <-------
rm(list=ls()[!ls() %in% c("DC_theme_generator")])

DF <- data.frame(expr=NA,time=NA,group=NA)
DF <- DF[-1,]
testname <- "bmark_base_V_stringi_"
for (f in list.files(path = "bmark/",pattern = testname)){
  tmp <- readRDS(paste0("bmark/",f))
  tmp$group <- rep(f,length(tmp[,1]))
  DF <- rbind(DF,tmp)
}
DF$group <- gsub(testname, "", DF$group, fixed = T)
DF$group <- gsub(".RDS", "", DF$group, fixed = T)
rm(tmp,testname,f)
DF$time <- DF$time/1e+06
DF$group[DF$group == "296602"] <- "300000"
DF <- DF[DF$group != "100", ]

customtheme <- DC_theme_generator(type = 'L',
                                  legend = 'F',
                                  ticks = 'out',
                                  x.axis.angle = 45,
                                  hjust = 1,
                                  vjust = 1,
                                  fontsize.cex = 1.2,
                                  ax.fontstyle = "italic")

ggplot(data = DF, 
       aes(x = expr, 
           y = time, 
           group = interaction(group,
                               expr)))+
  stat_summary(aes(col=expr),
               fun = median,
               geom = "crossbar", 
               width = 0.75,
               size = 0.25)+
  geom_dotplot(binaxis = "y",
               stackdir = "center", 
               position = "dodge",
               aes(fill=expr), 
               dotsize = 0.5,
               stroke=0.5)+
  customtheme+
  ylab("Time (µs)")+
  xlab("Method")+
  facet_wrap(~group ,
             nrow = 1,
             drop=T,
             scales = "free")+
  ggtitle("String operations")+
  expand_limits(y = 0)+
  scale_y_continuous(n.breaks = 4)

ggsave("baseVstringi.pdf", height = 4, width = 4)

#-------> data.frame vs data.table <-------
rm(list=ls()[!ls() %in% c("DC_theme_generator")])

DF <- data.frame(expr=NA,time=NA,group=NA)
DF <- DF[-1,]
testname <- "bmark_data.frame_vs_data.table_"
for (f in list.files(path = "bmark/",pattern = testname)){
  tmp <- readRDS(paste0("bmark/",f))
  tmp$group <- rep(f,length(tmp[,1]))
  DF <- rbind(DF,tmp)
}
DF$group <- gsub(testname, "", DF$group, fixed = T)
DF$group <- gsub(".RDS", "", DF$group, fixed = T)
rm(tmp,testname,f)

DF$time <- DF$time/1e+06

DF <- DF[DF$expr != "data.tablex4", ]
DF$expr <- gsub("data.tablex1","Data Table",DF$expr)
DF$expr <- gsub("data.frame","Data Frame",DF$expr)
DF <- DF[DF$group != "10000", ]


customtheme <- DC_theme_generator(type = 'L',
                                  legend = 'F',
                                  ticks = 'out',
                                  x.axis.angle = 45,
                                  hjust = 1,
                                  vjust = 1,
                                  fontsize.cex = 1.2,
                                  ax.fontstyle = "italic")

ggplot(data = DF, 
       aes(x = expr, 
           y = time, 
           group = interaction(group,
                               expr)))+
  stat_summary(aes(col=expr),
               fun = median,
               geom = "crossbar", 
               width = 0.75,
               size = 0.25)+
  geom_dotplot(binaxis = "y",
               stackdir = "center", 
               position = "dodge",
               aes(fill=expr), 
               dotsize = 0.5,
               stroke=0.5)+
  customtheme+
  ylab("Time (µs)")+
  xlab("Method")+
  facet_wrap(~group ,
             nrow = 1,
             drop=T,
             scales = "free")+
  ggtitle("Using function paramenters")+
  expand_limits(y = 0)+
  scale_y_continuous(n.breaks = 4)

ggsave("DF vs DT.pdf",
       height = 4,
       width = 4)


#-------> readRDS <-------
rm(list=ls()[!ls() %in% c("DC_theme_generator")])
results <- readRDS("results/results_parallel_readRDS.RDS")

DF <- readRDS("bmark/bmark_parallel_readRDS.RDS")
DF$time <- DF$time/1e+06

customtheme <- DC_theme_generator(type = 'L',
                                  legend = 'F',
                                  ticks = 'out',
                                  x.axis.angle = 45,
                                  hjust = 1,
                                  vjust = 1,
                                  fontsize.cex = 1.2,
                                  ax.fontstyle = "italic")

ggplot(data = DF, 
       aes(x = expr, 
           y = time, 
           ))+
  stat_summary(aes(col=expr),
               fun = median,
               geom = "crossbar", 
               width = 0.75,
               size = 0.25)+
  geom_dotplot(binaxis = "y",
               stackdir = "center", 
               position = "dodge",
               aes(fill=expr), 
               dotsize = 0.5,
               stroke=0.5)+
  customtheme+
  ylab("Time (µs)")+
  xlab("Method")+
  ggtitle("Saving RDS")+
  expand_limits(y = 0)+
  scale_y_continuous(n.breaks = 4)

ggsave("readRDS.pdf",
       height = 4,
       width = 2.5)

#-------> saveRDS <-------
rm(list=ls()[!ls() %in% c("DC_theme_generator")])

DF <- data.frame(expr=NA,time=NA,group=NA)
DF <- DF[-1,]
testname <- "bmark_parallel_saveRDS_"
for (f in list.files(path = "bmark/",pattern = testname)){
  tmp <- readRDS(paste0("bmark/",f))
  tmp$group <- rep(f,length(tmp[,1]))
  DF <- rbind(DF,tmp)
}

DF$time <- DF$time/1e+06
DF$group <- gsub(testname, "", DF$group, fixed = T)
DF$group <- gsub(".RDS", "", DF$group, fixed = T)
DF <- DF[DF$group != "10000",]

rm(tmp,testname,f)

customtheme <- DC_theme_generator(type = 'L',
                                  legend = 'F',
                                  ticks = 'out',
                                  x.axis.angle = 45,
                                  hjust = 1,
                                  vjust = 1,
                                  fontsize.cex = 1.2,
                                  ax.fontstyle = "italic")

ggplot(data = DF, 
       aes(x = expr, 
           y = time, 
           group = interaction(group,
                               expr)))+
  stat_summary(aes(col=expr),
               fun = median,
               geom = "crossbar", 
               width = 0.75,
               size = 0.25)+
  geom_dotplot(binaxis = "y",
               stackdir = "center", 
               position = "dodge",
               aes(fill=expr), 
               dotsize = 0.5,
               stroke=0.5)+
  customtheme+
  ylab("Time (µs)")+
  xlab("Method")+
  facet_wrap(~group ,
             nrow = 1,
             drop=T,
             scales = "free")+
  ggtitle("Using function paramenters")+
  expand_limits(y = 0)+
  scale_y_continuous(n.breaks = 4)

ggsave("saveRDS.pdf",
       height = 4,
       width = 4)


#-------> reading TSV <-------
rm(list=ls()[!ls() %in% c("DC_theme_generator")])

DF <- data.frame(expr=NA,time=NA,group=NA)
DF <- DF[-1,]
testname <- "bmark_readingTSV_"
for (f in list.files(path = "bmark/",pattern = testname)){
  tmp <- readRDS(paste0("bmark/",f))
  tmp$group <- rep(f,length(tmp[,1]))
  DF <- rbind(DF,tmp)
}

DF$time <- DF$time/1e+06
DF$group <- gsub(testname, "", DF$group, fixed = T)
DF$group <- gsub("_d.tsv.RDS", "", DF$group, fixed = T)
DF$group <- gsub("01", "00", DF$group)
DF <- DF[DF$group != "100",]

rm(tmp,testname,f)

customtheme <- DC_theme_generator(type = 'L',
                                  legend = 'F',
                                  ticks = 'out',
                                  x.axis.angle = 45,
                                  hjust = 1,
                                  vjust = 1,
                                  fontsize.cex = 1.2,
                                  ax.fontstyle = "italic")

ggplot(data = DF, 
       aes(x = expr, 
           y = time, 
           group = interaction(group,
                               expr)))+
  stat_summary(aes(col=expr),
               fun = median,
               geom = "crossbar", 
               width = 0.75,
               size = 0.25)+
  geom_dotplot(binaxis = "y",
               stackdir = "center", 
               position = "dodge",
               aes(fill=expr), 
               dotsize = 0.5,
               stroke=0.5)+
  customtheme+
  ylab("Time (µs)")+
  xlab("Method")+
  facet_wrap(~group ,
             nrow = 1,
             drop=T,
             scales = "free")+
  ggtitle("Using function paramenters")+
  expand_limits(y = 0)+
  scale_y_continuous(n.breaks = 4)

ggsave("readingTSV.pdf",
       height = 4,
       width = 5)


#-------> fixed and unlist <-------
rm(list=ls()[!ls() %in% c("DC_theme_generator")])

DF <- data.frame(expr=NA,time=NA,group=NA)
DF <- DF[-1,]
testname <- "bmark_fixed_and_unlist_"
for (f in list.files(path = "bmark/",pattern = testname)){
  tmp <- readRDS(paste0("bmark/",f))
  tmp$group <- rep(f,length(tmp[,1]))
  DF <- rbind(DF,tmp)
}

DF$time <- DF$time/1e+06
DF$group <- gsub(testname, "", DF$group, fixed = T)
DF$group <- gsub(".RDS", "", DF$group, fixed = T)
DF$group[DF$group == "1e+05"] <- "100000"
rm(tmp,testname,f)

customtheme <- DC_theme_generator(type = 'L',
                                  legend = 'F',
                                  ticks = 'out',
                                  x.axis.angle = 45,
                                  hjust = 1,
                                  vjust = 1,
                                  fontsize.cex = 1.2,
                                  ax.fontstyle = "italic")

ggplot(data = DF, 
       aes(x = expr, 
           y = time, 
           group = interaction(group,
                               expr)))+
  stat_summary(aes(col=expr),
               fun = median,
               geom = "crossbar", 
               width = 0.75,
               size = 0.25)+
  geom_dotplot(binaxis = "y",
               stackdir = "center", 
               position = "dodge",
               aes(fill=expr), 
               dotsize = 0.5,
               stroke=0.5)+
  customtheme+
  ylab("Time (µs)")+
  xlab("Method")+
  facet_wrap(~group ,
             nrow = 1,
             drop=T,
             scales = "free")+
  ggtitle("Using function paramenters")+
  expand_limits(y = 0)+
  scale_y_continuous(n.breaks = 4)

ggsave("fixed_and_unlist.pdf",
       height = 4,
       width = 4)
