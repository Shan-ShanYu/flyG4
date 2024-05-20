rm(list = ls());gc();rm(list = ls())
Num = "003.7."
setwd("/home/yuss/flyG4/result/data.reliability")
kc_top1000 <- fread("003.7.KC_top1000_summits.bed") %>% as.data.frame()
s2_top1000 <- fread("003.7.S2_top1000_summits.bed") %>% as.data.frame()

kc_top1000 <- kc_top1000[,1:2]
colnames(kc_top1000) <- c("chr","peak")
chr.list <- c("2L","2R","3L","3R","4","X","Y")
library(tidyverse)
kc_top1000_filter <- filter(kc_top1000, chr %in% chr.list)
kc_top1000_filter$start <- kc_top1000_filter$peak-100
kc_top1000_filter$end <- kc_top1000_filter$peak+100
kc_top1000_filter <- kc_top1000_filter[,-2]

s2_top1000 <- s2_top1000[,1:2]
colnames(s2_top1000) <- c("chr","peak")
s2_top1000_filter <- filter(s2_top1000, chr %in% chr.list)
s2_top1000_filter$start <- s2_top1000_filter$peak-100
s2_top1000_filter$end <- s2_top1000_filter$peak+100
s2_top1000_filter <- s2_top1000_filter[,-2]

write.table(kc_top1000_filter,file = paste0(Num,"KC_top1000_summits100_filter.bed"),
            sep = '\t',col.names = F,row.names = F,quote = F)
write.table(s2_top1000_filter,file = paste0(Num,"S2_top1000_summits100_filter.bed"),
            sep = '\t',col.names = F,row.names = F,quote = F)

kc_top1000_narrowPeak <-fread("003.7.KC_top1000.narrowPeak") %>%as.data.frame()
kc_top1000_narrowPeak <- kc_top1000_narrowPeak[,1:3]
colnames(kc_top1000_narrowPeak) <- c("chr","start","end")
kc_top1000_narrowPeak_filter <- filter(kc_top1000_narrowPeak, chr %in% chr.list)
s2_top1000_narrowPeak <- fread("003.7.S2_top1000.narrowPeak") %>% as.data.frame()
s2_top1000_narrowPeak <- s2_top1000_narrowPeak[,1:3]
colnames(s2_top1000_narrowPeak) <- c("chr","start","end")
s2_top1000_narrowPeak_filter <- filter(s2_top1000_narrowPeak, chr %in% chr.list)

write.table(kc_top1000_narrowPeak_filter,file = paste0(Num,"KC_top1000_filter.narrowPeak"),
            sep = '\t',col.names = F,row.names = F,quote = F)
write.table(s2_top1000_narrowPeak_filter,file = paste0(Num,"S2_top1000_filter.narrowPeak"),
            sep = '\t',col.names = F,row.names = F,quote = F)
# 003.7.S2_top1000_filter.narrowPeak.fa 958 sequences, 628748 residues
# TKYTKYTKYTGYTKYTGYTKY 3851
# CACACACACACACACACRCAC	3774
# AGMGAGAGMGAGMGAGAGMGV 3050
# VMAAMDNCAAMRRCARCRACA	2799
# BRTRYRTWCRTAYATATATRTATDTATAT	682
# BCGCYGCCGVCGYCGMYGCCG	1503 
# 003.7.KC_top1000_filter.narrowPeak.fa 931 sequences, 592245 residues
# CACACACACACACACACACAC	3291
S2.density <- data.frame(type=c("motif 1","motif 2","motif 3","motif 4","motif 5","motif 6"),
                              occurences=c(3851,3774,3050,2799,682,1503),
                              residues=c(rep(628748,6)),
                              GC=c("more","more","more","more","less","less"))
S2.density$density <- round(S2.density$occurences/S2.density$residues,5)
S2.density$type <- factor(S2.density$type, levels = c("motif 1","motif 2","motif 3","motif 4","motif 5","motif 6"))
ggplot(data = S2.density,aes(x=type,y=density,fill=GC))+
  geom_bar(stat="identity",width=0.7,position='dodge')+
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("black", "#F9A500")) +
  ylab("Density (occurence/residue)") +
  cowplot::theme_half_open() +
  theme(axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        legend.position = "none") +
  coord_cartesian(ylim = c(0,0.007))

# 003.7.KC_top1000_filter.narrowPeak.fa 931 sequences, 592245 residues
# CACACACACACACACACACAC	3291
# RMARCARMARMARCARCARMA	3334
# CTCTCKCTCTCKCTCTCTCKCTCTCTY	3989
# MAAAANVRMAAMAAMAACAAA	1304
# TRTRTGTRYRTRTGTRTGTGYRTDYRT 3135
# NGNGVGWGMGAGCGAGASRG 2453
# CGSCGGCRGCGRCGNCGVCRRCGVMGVCG 2544
# GCAGCASCANCASCA 1758
# GYRYSTGTGTGTGTG 2632 
KC.density <- data.frame(type=c("motif 1","motif 2","motif 3","motif 4","motif 5","motif 6","motif 7","mofit 8","motif 9"),
                              occurences=c(3291,3334,3989,1304,3135,2453,2544,1758,2632),
                              residues=c(rep(592245,9)),
                              GC=c("more","more","more","less","more","more","more","less","more"))
KC.density$density <- round(KC.density$occurences/KC.density$residues,5)
KC.density$type <- factor(KC.density$type, levels = c("motif 1","motif 2","motif 3","motif 4","motif 5","motif 6","motif 7","mofit 8","motif 9"))
ggplot(data = KC.density,aes(x=type,y=density,fill=GC))+
  geom_bar(stat="identity",width=0.7,position='dodge')+
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("black", "#F9A500")) +
  ylab("Density (occurence/residue)") +
  cowplot::theme_half_open() +
  theme(axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        legend.position = "none") +
  coord_cartesian(ylim = c(0,0.007))
