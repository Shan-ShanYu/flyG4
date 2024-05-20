rm(list = ls());gc();rm(list = ls())
Num = "003.3."
setwd("/home/yuss/flyG4/result/data.reliability")
library(tidyverse)
kc.rep1 <- fread("KC-G4-rep1_R1_summits.bed") %>% as.data.frame()
kc.rep1 <- kc.rep1[,1:2]
colnames(kc.rep1) <- c("chr","peak")
chr.list <- c("2L","2R","3L","3R","4","X","Y")
kc.rep1.filter <- filter(kc.rep1, chr %in% chr.list)
kc.rep1.filter$start <- kc.rep1.filter$peak-100
kc.rep1.filter$end <- kc.rep1.filter$peak+100
kc.rep1.filter <- kc.rep1.filter[,-2]

kc.rep2 <- fread("KC-G4-rep2_R1_summits.bed") %>% as.data.frame()
kc.rep2 <- kc.rep2[,1:2]
colnames(kc.rep2) <- c("chr","peak")
kc.rep2.filter <- filter(kc.rep2, chr %in% chr.list)
kc.rep2.filter$start <- kc.rep2.filter$peak-100
kc.rep2.filter$end <- kc.rep2.filter$peak+100
kc.rep2.filter <- kc.rep2.filter[,-2]

s2.rep1 <- fread("S2-G4-rep1_R1_summits.bed") %>% as.data.frame()
s2.rep1 <- s2.rep1[,1:2]
colnames(s2.rep1) <- c("chr","peak")
s2.rep1.filter <- filter(s2.rep1, chr %in% chr.list)
s2.rep1.filter$start <- s2.rep1.filter$peak-100
s2.rep1.filter$end <- s2.rep1.filter$peak+100
s2.rep1.filter <- s2.rep1.filter[,-2]

s2.rep2 <- fread("S2-G4-rep2_R1_summits.bed") %>% as.data.frame()
s2.rep2 <- s2.rep2[,1:2]
colnames(s2.rep2) <- c("chr","peak")
s2.rep2.filter <- filter(s2.rep2, chr %in% chr.list)
s2.rep2.filter$start <- s2.rep2.filter$peak-100
s2.rep2.filter$end <- s2.rep2.filter$peak+100
s2.rep2.filter <- s2.rep2.filter[,-2]

write.table(kc.rep1.filter, file = paste0(Num,"KC.rep1.summits100.filter.bed"),
            sep = '\t',col.names = F,row.names = F,quote = F)
write.table(kc.rep2.filter, file = paste0(Num,"KC.rep2.summits100.filter.bed"),
            sep = '\t',col.names = F,row.names = F,quote = F)
write.table(s2.rep1.filter, file = paste0(Num,"S2.rep1.summits100.filter.bed"),
            sep = '\t',col.names = F,row.names = F,quote = F)
write.table(s2.rep2.filter, file = paste0(Num,"S2.rep2.summits100.filter.bed"),
            sep = '\t',col.names = F,row.names = F,quote = F)
#KC.rep1 Database contains 46126 sequences, 9225200 residues
#1.ACACACACACACACA 17125 motif occurences 
#2.GCTGYTGCTGCTGCTGCTGCT 52648 motif occurences
#3.MAAMAVMAAWNVMAA 10774 motif occurences
#4.GYGGWGGWGGWGGHGGWGGDG 37424 motif occurences
#5.STCSTCGTCGTCGTC 8319 motif occurences
#6(7).GAGAGAGAGAGAGAGAGAGRGAGWGRGAG 14021 motif occurences
KC.rep1.density <- data.frame(type=c("motif 1","motif 2","motif 3","motif 4","motif 5","motif 6"),
                                 occurences=c(17125,52648,10774,37424,8319,14021),
                                 residues=c(rep(9225200,6)),
                              GC=c("more","more","less","more","less","more"))
KC.rep1.density$density <- round(KC.rep1.density$occurences/KC.rep1.density$residues,5)
library(ggplot2)
ggplot(data = KC.rep1.density,aes(x=type,y=density,fill=GC))+
  geom_bar(stat="identity",width=0.8,position='dodge')+ #,color = "black"边框是黑色
  scale_fill_manual(values=c("black", "#F9A500")) +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Density (occurence/residue)") +
  cowplot::theme_half_open() +
  theme(axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        legend.position = "none") + #删除图例
  coord_cartesian(ylim = c(0,0.006))
#KC.rep2 Database contains 50066 sequences, 10013200 residues 
#1.GCARCAGCAGCARCAGCAGCA 50633 motif occurences
#2.MMHCACACACACACACACACA 21826 motif occurences
#3.CTCCTCCTCCTCCTC 19599 motif occurences
#4.MAARRMMAAMAAMAA 6961 motif occurences
#5.GSHGVYGGHGRTGGHGRTGGH 40250 motif occurences
KC.rep2.density <- data.frame(type=c("motif 1","motif 2","motif 3","motif 4","motif 5"),
                              occurences=c(50633,21826,19599,6961,40250),
                              residues=c(rep(9225200,5)),
                              GC=c("more","more","more","less","more"))
KC.rep2.density$density <- round(KC.rep2.density$occurences/KC.rep2.density$residues,5)
library(ggplot2)
ggplot(data = KC.rep2.density,aes(x=type,y=density,fill=GC))+
  geom_bar(stat="identity",width=0.8,position='dodge')+
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values=c("black", "#F9A500")) +
  ylab("Density (occurence/residue)") +
  cowplot::theme_half_open() +
  theme(axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        legend.position = "none") +
  coord_cartesian(ylim = c(0,0.006))
#S2.rep1 Database contains 41703 sequences, 8340600 residues
#1.CAGCAGCAGCAGCA	28557 motif occurences
#2.TGYGTGYGTGTGTGTGTGTGTGTGTGTGT 23952 motif occurences
#3.AAAAMAAMAARARAMAAAARR	5414 motif occurences
#4.GKKGGHGGTGGDGGTGGNBGKGGBGKKGG 36074 motif occurences
#5.SYTGGCCA 9769 motif occurences
#6.GCRGCAACARCARCR 19480 motif occurences
#7.TGGCTGGCTGGMTGGCTGRCT 6292 motif occurences
#8.CABATATAYRTATGTATG 2172 motif occurences
#9.GATGVDGMTGMDGATGVDGAD 16919 motif occurences
S2.rep1.density <- data.frame(type=c("motif 1","motif 2","motif 3","motif 4","motif 5","motif 6","motif 7","mofit 8","motif 9"),
                              occurences=c(28557,23952,5414,36074,9769,19480,6292,2172,16919),
                              residues=c(rep(8340600,9)),
                              GC=c("more","more","less","more","less","more","more","less","more"))
S2.rep1.density$density <- round(S2.rep1.density$occurences/S2.rep1.density$residues,5)
S2.rep1.density$type <- factor(S2.rep1.density$type, levels = c("motif 1","motif 2","motif 3","motif 4","motif 5","motif 6","motif 7","mofit 8","motif 9"))
ggplot(data = S2.rep1.density,aes(x=type,y=density,fill=GC))+
  geom_bar(stat="identity",width=0.8,position='dodge')+
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("black", "#F9A500")) +
  ylab("Density (occurence/residue)") +
  cowplot::theme_half_open() +
  theme(axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        legend.position = "none") +
  coord_cartesian(ylim = c(0,0.005))
#S2.rep2 Database contains 46116 sequences, 9223200 residues
#1.TGYTGCTGCTGCTGCTGCTGC 44724 motif occurences
#2.ACACACACACACACACACACA	24274 motif occurences
#3.GAGGWGGAKGAGGAKGWGGAK 24404 motif occurences
#4.VRMMAAAAAAMAVAAAAHVRAAADVAAA	11586 motif occurences
#5.CTGCWGCTSCT 19427 motif occurences
#6.AWTGCAATTGSCAWT 4143 motif occurences
#7.TYYYCCTCTCTCYCTCTCTCTCTCK	13170 motif occurences
S2.rep2.density <- data.frame(type=c("motif 1","motif 2","motif 3","motif 4","motif 5","motif 6","motif 7"),
                              occurences=c(44724,24274,24404,11586,19427,4143,13170),
                              residues=c(rep(9223200,7)),
                              GC=c("more","more","more","less","more","less","more"))
S2.rep2.density$density <- round(S2.rep2.density$occurences/S2.rep2.density$residues,5)
S2.rep2.density$type <- factor(S2.rep2.density$type, levels = c("motif 1","motif 2","motif 3","motif 4","motif 5","motif 6","motif 7"))
ggplot(data = S2.rep2.density,aes(x=type,y=density,fill=GC))+
  geom_bar(stat="identity",width=0.8,position='dodge')+
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("black","#F9A500")) +
  ylab("Density (occurence/residue)") +
  cowplot::theme_half_open() +
  theme(axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        legend.position = "none") +
  coord_cartesian(ylim = c(0,0.005))

