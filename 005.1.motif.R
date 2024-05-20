####用seqLogo可视乎motif####
#PPM-PFM-PWM
# BiocManager::install("seqLogo")
library(seqLogo)
#1.输入motif对应的PPM矩阵PPM（Position Probability Matrix），也称为位置概率矩阵
data <- read.table("~/flyG4/1.txt.txt",header = F,sep = "\t",row.names = 1)
#2.计算PFM
##PPM矩阵就是将PFM矩阵中的频数转化成频率，除以每列的总和就可以了Position Frequency Matrix（PFM）。
ppm <- sapply(1:ncol(data), function(t){data[[t]]/sum(data[[t]])})
#3.位置权重矩阵position weight matrix（PWM）
p <- makePWM(ppm)
seqLogo(p)　　　　　　

####用ggseqlogo修改配色####
#*KC motif---------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())#清空
Num = "005.1."
# BiocManager::install("ggseqlogo")
library(ggseqlogo)
# install.packages("ggplot2")
library(ggplot2)
library(magrittr)

setwd("/home/yuss/flyG4/result/TopMotif/")
filename <- list.files("/home/yuss/flyG4/result/TopMotif/",pattern = "^kc")
a <- read.table("/home/yuss/flyG4/result/TopMotif/kc.FBgn0001325_3.txt",header = F,sep = "\t") %>% as.matrix()
# 创建一个空列表来存储矩阵
matrix_list <- list()
for (i in 1:length(filename)){
  var_name <- gsub('.txt', '',filename[i])
  data_matrix <- as.matrix(read.table(filename[i], sep = '\t', header = F)) ##assign()函数将一个读取的数据框对象分配给先前定义的变量名 var_name
  rownames(data_matrix) <- c("A", "C", "G", "T")
  matrix_list[[var_name]] <- data_matrix
  }
ggseqlogo(matrix_list,method="bits",col_scheme="nucleotide",facet = "wrap",ncol = 4)

# 翻转
flip <- matrix_list[["kc.FBgn0003870_4"]]
flip_matrix <- flip[nrow(flip):1,ncol(flip):1]
rownames(flip_matrix) <- c("A", "C", "G", "T")
matrix_list[["kc.FBgn0003870_4"]] <- flip_matrix

flip <- matrix_list[["kc.FBgn0005630_14"]]
flip_matrix <- flip[nrow(flip):1,ncol(flip):1]
rownames(flip_matrix) <- c("A", "C", "G", "T")
matrix_list[["kc.FBgn0005630_14"]] <- flip_matrix
names(matrix_list) <- c("Motif 1","Motif 2","Motif 3","Motif 4","Motif 5","Motif 6","Motif 7","Mofit 8")
ggseqlogo(matrix_list,method="bits",col_scheme="nucleotide",facet = "wrap",ncol = 4)
ggsave(filename = paste0("/home/yuss/flyG4/result/TopMotif/Picture/",Num,"KCMotif.pdf"),
       device = "pdf",width = 10,height = 3)

#*KC density---------------------------------------------------------------------
# Database contains 37127 sequences, 3750405 residues
# OTF0397.1	Q8MR37_DROME_B1H	AWGGGCGTGGC 9.6e-62 2095 motif occurences
# OTF0352.1	Q7K9G4_DROME_B1H	BGYGGGGGGKS 9.6e-54 2629 motif occurences
# FBgn0003870_4	ttk-PA_SANGER_5	CMACCCCTA 7.1e-35 530 motif occurences  (相反链)
# FBgn0013469	klu_SANGER_10	TGYGKGGGTGK 5.4e-32 3280 motif occurences
# FBgn0005630_15	lola-PL_SANGER_2.5	KGTGGGKCA 5.1e-19 1486 motif occurences
# FBgn0001325_3	Kr_SANGER_5	RARGGGTW 1.1e-9 400 motif occurences
# FBgn0005630_14	lola-PF_SANGER_5	MAACTCCAYY 6.1e-6 (相反链) 1445 motif occurences 
# GGANGNGGAKGHGGA	MEME-3	GGANGNGGAKGHGGA 2.8e-4 8070 motif occurences

kc.density <- data.frame(name=c("FBgn0001325_3","FBgn0003870_4","FBgn0005630_14","FBgn0005630_15",
                                "FBgn0013469","GGANGNGGAKGHGGA","OTF0352.1","OTF0397.1"),
                         type=c("Motif 1","Motif 2","Motif 3","Motif 4","Motif 5","Motif 6","Motif 7","Mofit 8"),
                         occurences=c(400,530,1445,1486,3280,8070,2629,2095),
                         residues=c(rep(3750405,8)))
kc.density$density <- round(kc.density$occurences/kc.density$residues,5)
kc.density$type <- factor(kc.density$type, levels = c("Motif 1","Motif 2","Motif 3","Motif 4","Motif 5","Motif 6","Motif 7","Mofit 8"))
ggplot(data = kc.density,aes(x=type,y=density))+
  geom_bar(stat="identity",width=0.7,position='dodge',color = "black", fill = "#fdb863")+
  scale_y_continuous(expand = c(0,0)) +
  ylab("Density (occurence/residue)") +
  cowplot::theme_half_open() +
  theme(axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        legend.position = "none") +
  coord_cartesian(ylim = c(0,0.003))
ggsave(filename = paste0("/home/yuss/flyG4/result/TopMotif/Picture/",Num,"KCMotif.Density.pdf"),
       device = "pdf",width = 6,height = 3.2)

#*S2---------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())#清空
Num = "005.1."

library(ggseqlogo)
library(ggplot2)
library(magrittr)
setwd("/home/yuss/flyG4/result/TopMotif/")
filename <- list.files("/home/yuss/flyG4/result/TopMotif/",pattern = "^s2")

# 创建一个空列表来存储矩阵
matrix_list <- list()
for (i in 1:length(filename)){
  var_name <- gsub('.txt', '',filename[i])
  data_matrix <- as.matrix(read.table(filename[i], sep = '\t', header = F)) ##assign()函数将一个读取的数据框对象分配给先前定义的变量名 var_name
  rownames(data_matrix) <- c("A", "C", "G", "T")
  matrix_list[[var_name]] <- data_matrix
}
ggseqlogo(matrix_list,method="bits",col_scheme="nucleotide",facet = "wrap",ncol = 4)
##把含有C的motif的矩阵行和列进行翻转
flip <- matrix_list[["s2.Ci_SANGER_5"]]
flip_matrix <- flip[nrow(flip):1,ncol(flip):1]
rownames(flip_matrix) <- c("A", "C", "G", "T")
matrix_list[["s2.Ci_SANGER_5"]] <- flip_matrix

flip <- matrix_list[["s2.OTF0267.1"]]
flip_matrix <- flip[nrow(flip):1,ncol(flip):1]
rownames(flip_matrix) <- c("A", "C", "G", "T")
matrix_list[["s2.OTF0267.1"]] <- flip_matrix

flip <- matrix_list[["s2.ttk-PA_SANGER_5"]]
flip_matrix <- flip[nrow(flip):1,ncol(flip):1]
rownames(flip_matrix) <- c("A", "C", "G", "T")
matrix_list[["s2.ttk-PA_SANGER_5"]] <- flip_matrix
names(matrix_list) <- c("Motif 1","Motif 2","Motif 3","Motif 4","Motif 5","Motif 6","Motif 7","Mofit 8")
ggseqlogo(matrix_list,method="bits",col_scheme="nucleotide",facet = "wrap",ncol = 4)
ggsave(filename = paste0("/home/yuss/flyG4/result/TopMotif/Picture/",Num,"S2Motif.pdf"),
       device = "pdf",width = 10,height = 3)
# Database contains 29667 sequences, 2996805 residues
# FBgn0086910_2	l(3)neo38_SOLEXA_2.5	DGKGGGKGGGGGDGD	1.2e-49 5028 motif occurences
##### FBgn0000568	Eip75B_SANGER_5	TAWDTRGGTCA	5.1e-6 
# FBgn0003870_4	ttk-PA_SANGER_5	CMACCCCTA	7.7e-27 435 motif occurences
# FBgn0004859	Ci_SANGER_5	RGACCACCCAC	2.2e-11 1398 motif occurences
# OTF0013.1	A4IJ80_DROME_B1H	AGTGGGCGKGR	1.5e-4 2029 motif occurences
# FBgn0020309_2	crol-F7-16_SOLEXA	GGGBGVGGGGGGGDD	3.1e-12 4669 motif occurences
# OTF0267.1	OPA_DROME_B1H	GMCCCCCCGCT	8.7e-12 1196 motif occurences 
# OTF0470.1	Q9VN10_DROME_B1H	KGGGCGTGA	4.3e-5 1308 motif occurences 
# 28-AGGATGTGGA	STREME-28	AGGATGTGGA	8.8e-22 2122 motif occurences
## 27-AGGAGCTG	STREME-27	AGGAGCTG	8.2e-4 440 motif occurences
## 9-GCCCCAAAA 2.9e-002 1230 motif occurences

S2.density <- data.frame(name=c("28-AGGATGTGGA","A4IJ80_DROME_B1H","Ci_SANGER_5","crol-F7-16_SOLEXA",
                                "l(3)neo38_SOLEXA_2.5","OTF0267.1","OTF0470.1","ttk-PA_SANGER_5"),
                         type=c("Motif 1","Motif 2","Motif 3","Motif 4","Motif 5","Motif 6","Motif 7","Motif 8"),
                         occurences=c(2122,2029,1398,4669,5028,1196,1308,435),
                         residues=c(rep(2996805,8)))
S2.density$density <- round(S2.density$occurences/S2.density$residues,5)
S2.density$type <- factor(S2.density$type, levels = c("Motif 1","Motif 2","Motif 3","Motif 4","Motif 5","Motif 6","Motif 7","Motif 8"))
ggplot(data = S2.density,aes(x=type,y=density))+
  geom_bar(stat="identity",width=0.7,position='dodge',color = "black", fill = "#fdb863")+
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = "yello") +
  ylab("Density (occurence/residue)") +
  cowplot::theme_half_open() +
  theme(axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        legend.position = "none") +
  coord_cartesian(ylim = c(0,0.002))
ggsave(filename = paste0("/home/yuss/flyG4/result/TopMotif/Picture/",Num,"S2Motif.Density.pdf"),
       device = "pdf",width = 6,height = 3.2)
