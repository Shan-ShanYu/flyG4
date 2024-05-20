rm(list = ls());gc();rm(list = ls())#清空
Num = "001.4."
#install.packages("eulerr")
library(eulerr)
VennDiag <- euler(c("G4-CUTTag" = 2564,"G4-Seq(PDS)" = 41424,
                    "G4-CUTTag&G4-Seq(PDS)" = 13839))
plot(VennDiag, counts = TRUE, font=3, cex=1, alpha=1,quantities = TRUE,lwd =3.5,
          labels=c("G4-CUT&Tag","G4-Seq(PDS)"),
          label.col = "white",fill="white",
          col = c('#e48385','#4974a5'))

VennDiag <- euler(c("G4-CUTTag" = 3062,"G4-Seq(K+)" = 2246,"G4-Seq(PDS)" = 29090,
                    "G4-CUTTag&G4-Seq(K+)" = 27, "G4-CUTTag&G4-Seq(PDS)" = 9047, "G4-Seq(PDS)&G4-Seq(K+)" = 12860,
                    "G4-CUTTag&G4-Seq(K+)&G4-Seq(PDS)" = 4266))
plot(VennDiag, counts = TRUE, font=3, cex=1, alpha=1,quantities = TRUE,lwd =3.5,
     labels=c("G4-CUT&Tag","G4-Seq(K+)","G4-Seq(PDS)"),
     label.col = "white",fill="white",
     col = c('#e48385','#4974a5','#66c2a5'))

##peak四个合并，k+和PDS合并，再将两者取交集
VennDiag <- euler(c("G4-CUTTag" = 118095,"G4-Seq" = 7601,
                    "G4-CUTTag&G4-Seq" = 67061))
p <- plot(VennDiag, counts = TRUE, font=3, cex=1, alpha=1,quantities = TRUE,lwd =3.5,
     labels=c("G4-CUT&Tag","G4-Seq"),
     label.col = "white",fill="white",
     col = c('#e48385','#4974a5'))
ggsave(p,filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"PeakIntersectG4seq.pdf"),
       device = "pdf",width = 3.8,height = 3)

##peak两个重复之间取交集，再将kc s2的peak合并，k+和PDS合并，再将两者取交集
VennDiag <- euler(c("G4-CUTTag" = 80439-30072,"G4-Seq" = 74662-30072,
                    "G4-CUTTag&G4-Seq" = 30072))
p <- plot(VennDiag, counts = TRUE, font=3, cex=1, alpha=1,quantities = TRUE,lwd =3.5,
          labels=c("G4-CUT&Tag","G4-Seq"),
          label.col = "white",fill="white",
          col = c('#e48385','#4974a5'))
ggsave(p,filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"KCS2CatPeak.IntersectG4seq.pdf"),
       device = "pdf",width = 3.8,height = 3.8)

#### 每个细胞系的peak与PQS取交集 ####
在一个细胞系kc167中，当PQS与kc两个重复分别都有交集才认定为eG4，即两个重复必须有交集，且与PQS也交集，被认定为eG4，如果两个重复没有交集，被认定为non-eG4
# (base) yuss@ubuntu:~/flyG4/result/LucyCherbas.GR.2010.RNAseq$ bedtools intersect -a "/home/yuss/flyG4/result/PQS/001.2.KC-G4-rep2_R1_peaks.narrowPeak" -b "/home/yuss/flyG4/result/PQS/001.2.KC-G4-rep1_R1_peaks.narrowPeak" | grep -E "2L|2R|3L|3R" | wc -l
# 36315
# (base) yuss@ubuntu:~/flyG4/result/LucyCherbas.GR.2010.RNAseq$  bedtools intersect -a "/home/yuss/flyG4/result/PQS/001.2.S2-G4-rep2_R1_peaks.narrowPeak" -b "/home/yuss/flyG4/result/PQS/001.2.S2-G4-rep1_R1_peaks.narrowPeak" | grep -E "2L|2R|3L|3R" | wc -l
# 34550
##kc
VennDiag <- euler(c("G4-CUTTag" = 43192-13981,"PQS" = 43917-13981,
                    "G4-CUTTag&PQS" = 13981))
p <- plot(VennDiag, counts = TRUE, font=3, cex=1, alpha=1,quantities = TRUE,lwd =3.5,
     labels=c("Kc G4-CUT&Tag","PQS"),
          label.col = "white",fill="white",
          col = c('#D76364','#5B4E4A'))
ggsave(p,filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"KcIntersectPeak.IntersectPQS.pdf"),
       device = "pdf",width = 3.8,height = 3.8)
##s2
VennDiag <- euler(c("G4-CUTTag" = 38659-12004,"PQS" = 43917-12004,
                    "G4-CUTTag&PQS" = 12004))
p <- plot(VennDiag, counts = FALSE, font=3, cex=1, alpha=1,quantities = TRUE,lwd =3.5,
     labels=c("S2 G4-CUT&Tag","PQS"),
     label.col = "white",fill="white",
     col = c('#5F97D3','#5B4E4A'))
ggsave(p, filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"S2IntersectPeak.IntersectPQS.pdf"),
       device = "pdf",width = 3.8,height = 3.8)

# bedtools intersect -a "/home/yuss/flyG4/result/PQS/001.2.KC-G4-rep2_R1_peaks.narrowPeak" -b "/home/yuss/flyG4/result/PQS/001.2.KC-G4-rep1_R1_peaks.narrowPeak" -wa | grep -E -w "2L|2R|3L|3R|^4|X|Y" | wc -l
# 43192
# bedtools intersect -a "/home/yuss/flyG4/result/PQS/001.2.S2-G4-rep2_R1_peaks.narrowPeak" -b "/home/yuss/flyG4/result/PQS/001.2.S2-G4-rep1_R1_peaks.narrowPeak" | grep -E -w "2L|2R|3L|3R|4|X|Y" | wc -l
# 38659

# bedtools intersect -a "/home/yuss/flyG4/result/PQS/001.2.KC-G4-rep1_R1_peaks.narrowPeak" -b "/home/yuss/flyG4/result/PQS/001.2.KC-G4-rep2_R1_peaks.narrowPeak" | grep -E -w "2L|2R|3L|3R|^4|X|Y" > 001.4.KC-G4-repintersect.narrowPeak
# wc -l 001.4.KC-G4-repintersect.narrowPeak
# 43192 
# bedtools intersect -a 001.1.dmel.pqs.bed -b 001.4.KC-G4-repintersect.narrowPeak | wc -l
# 13981
# bedtools intersect -a "/home/yuss/flyG4/result/PQS/001.2.S2-G4-rep2_R1_peaks.narrowPeak" -b "/home/yuss/flyG4/result/PQS/001.2.S2-G4-rep1_R1_peaks.narrowPeak" | grep -E -w "2L|2R|3L|3R|^4|X|Y" > 001.4.S2-G4-repintersect.narrowPeak
# wc -l 001.4.S2-G4-repintersect.narrowPeak
# 38659 
# (base) yuss@ubuntu:~/flyG4/result/PQS$ bedtools intersect -a 001.1.dmel.pqs.bed -b 001.4.S2-G4-repintersect.narrowPeak | wc -l
# 12004
# (base) yuss@ubuntu:~/flyG4/result/PQS$ bedtools intersect -a 001.1.dmel.pqs.bed -b 001.4.KC-G4-repintersect.narrowPeak > 001.4.PqsIntersectKc.bed
# (base) yuss@ubuntu:~/flyG4/result/PQS$ bedtools intersect -a 001.4.PqsIntersectKc.bed -b 001.4.S2-G4-repintersect.narrowPeak | wc -l
# 9583
# (base) yuss@ubuntu:~/flyG4/result/PQS$ bedtools intersect -a 001.4.KC-G4-repintersect.narrowPeak -b 001.4.S2-G4-repintersect.narrowPeak | wc -l
# 32942

VennDiag <- euler(c("Kc167" = 5852,"S2" = 3296,"PQS" = 27515,
                    "Kc167&PQS" = 4398, "Kc167&S2" = 23359, "PQS&S2" = 2421,
                    "Kc167&PQS&S2" = 9583))
p <- plot(VennDiag, counts = TRUE, font=3, cex=1, alpha=1,quantities = TRUE,lwd =3.5,
     # labels=c("Kc167 G4-CUT&Tag","S2 G4-CUT&Tag","PQS"),
     label.col = "white",fill="white",
     col = c('#e48385','#4974a5','#66c2a5'))
ggsave(p, filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"KcS2PQSIntersect.pdf"),
       device = "pdf",width = 3.8,height = 3.8)

rm(list = ls());gc();rm(list = ls())#清空
Num = "001.4."
VennDiag <- euler(c("kcrep1" = 46461-43192,"kcrep2" = 50422-43192,
                    "kcrep1&kcrep2" = 43192))
p <- plot(VennDiag, counts = FALSE, font=3, cex=1, alpha=1,quantities = TRUE,lwd =3.5,
          labels=c("Kc167 Rep1","Kc167 Rep2"),
          label.col = "white",fill="white",
          col = c('#D76364','#D76364'))
ggsave(p, filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"KcIntersectPeak.pdf"),
       device = "pdf",width = 3.8,height = 3.8)

VennDiag <- euler(c("s2rep1" = 41922-38659,"s2rep2" = 46351-38659,
                    "s2rep1&s2rep2" = 38659))
p <- plot(VennDiag, counts = FALSE, font=3, cex=1, alpha=1,quantities = TRUE,lwd =3.5,
          labels=c("S2 Rep1","S2 Rep2"),
          label.col = "white",fill="white",
          col = c('#5F97D3','#5F97D3'))
ggsave(p, filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"S2IntersectPeak.pdf"),
       device = "pdf",width = 3.8,height = 3.8)
