#### DMSO ####
rm(list = ls());gc();rm(list = ls())
Num = "010.7."
readcounts <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.counts.txt") %>% as.data.frame()
rownames(readcounts) <- readcounts$geneid
con <- readcounts[,c(1:3,9:11)]
sample_name <- factor(colnames(con))
metadata <- data.frame(sample_name)
metadata$cell <- as.factor(rep(c("kc","s2"),each=3))

dds <- DESeqDataSetFromMatrix(countData = con,
                              colData = metadata,
                              design = ~ cell)
##2.dds标准化
dds <- DESeq(dds)
##3.获取标准化后的数据
normalized_counts <- counts(dds,normalized=T)
##4.有三个重复取平均
normalized_mean_counts = t(apply(normalized_counts, 1, function(a){tapply(a, metadata$cell, mean)}))

gene.kc.s2 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.2.gene.kc.s2.txt") %>% as.data.frame()
gene.kc.s2$kc.eG4 <- ifelse(gene.kc.s2$kc==1,"eG4","no eG4")
gene.kc.s2$s2.eG4 <- ifelse(gene.kc.s2$s2==2,"eG4","no eG4")
gene.kc.s2 <- subset(gene.kc.s2,chr=='2L'|chr=='2R'|chr=='3L'|chr=='3R'|chr=='X')
gene.kc.s2$kccon <- normalized_mean_counts[match(gene.kc.s2$id,rownames(normalized_mean_counts)),1]
gene.kc.s2$s2con <- normalized_mean_counts[match(gene.kc.s2$id,rownames(normalized_mean_counts)),2]

df <- na.omit(gene.kc.s2)
df$chr1 <- ifelse(df$chr=="X","X","A")
df$ratio <- df$kccon/df$s2con
df$ratio1 <- df$s2con/df$kccon
library(ggpubr)
my_comparisons = list(c("eG4","no eG4"))
ggplot(data = df,aes(x=kc.eG4,y=log2(ratio+1),fill=kc.eG4)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(4,4.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  # scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Control log2(Kc/S2)") +
  coord_cartesian(ylim = c(0,5)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title="Kc167 cell")

ggplot(data = df,aes(x=kc.eG4,y=log2(ratio+1),fill=kc.eG4)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  facet_grid(~chr1) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(4,4.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  # scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Control log2(Kc/S2+1)") +
  coord_cartesian(ylim = c(0,5)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title="Kc167 cell") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") 

non_finite_rows <- which(!is.finite(df$ratio))
# 如果需要，可以将非有限值替换为 NA
df$ratio[non_finite_rows] <- NA

# 删除包含缺失值的行
dfkc_clean <- df[complete.cases(df), ]

table(dfkc_clean$kc.eG4,dfkc_clean$chr1)
sample_data <- dfkc_clean[dfkc_clean$kc.eG4=="eG4"&dfkc_clean$chr1=="X",15]
wilcox.test(sample_data, mu = 1,alternative = "greater") #p-value < 2.2e-16
sample_data <- dfkc_clean[dfkc_clean$kc.eG4=="no eG4"&dfkc_clean$chr1=="X",15]
wilcox.test(sample_data, mu = 1,alternative = "greater") #p-value < 2.2e-16

df$ratio1 <- df$s2con/df$kccon
ggplot(data = df,aes(x=s2.eG4,y=log2(ratio1+1),fill=s2.eG4)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(4,4.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  # scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Control log2(S2/Kc+1)") +
  coord_cartesian(ylim = c(0,5)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title="S2 cell") 

ggplot(data = df,aes(x=s2.eG4,y=log2(ratio1+1),fill=s2.eG4)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  facet_grid(~chr1) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(4,4.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  # scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Control log2(S2/Kc+1)") +
  coord_cartesian(ylim = c(0,5)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title="S2 cell")+
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") 

non_finite_rows <- which(!is.finite(df$ratio1))
# 如果需要，可以将非有限值替换为 NA
df$ratio1[non_finite_rows] <- NA
# 删除包含缺失值的行
dfs2_clean <- df[complete.cases(df), ]

table(dfs2_clean$s2.eG4,dfs2_clean$chr1)
sample_data <- dfs2_clean[dfs2_clean$s2.eG4=="eG4"&dfs2_clean$chr1=="X",16]
wilcox.test(sample_data, mu = 1,alternative = "greater") #p-value 
sample_data <- dfs2_clean[dfs2_clean$s2.eG4=="no eG4"&dfs2_clean$chr1=="X",16]
wilcox.test(sample_data, mu = 1,alternative = "less") #p-value 

# 比较两组独立样本的中位数是否存在差异。
wilcox.test(dfkc_clean[dfkc_clean$s2.eG4=="eG4"&dfkc_clean$chr1=="X",13],dfkc_clean[dfkc_clean$kc.eG4=="eG4"&dfkc_clean$chr1=="X",12])
wilcox.test(df[df$s2.eG4=="eG4"&df$chr1=="X",13],df[df$kc.eG4=="eG4"&df$chr1=="X",12]) #p-value = 0.01416
wilcox.test(df[df$s2.eG4=="no eG4"&df$chr1=="X",13],df[df$kc.eG4=="no eG4"&df$chr1=="X",12]) #0.07191

#可视化 X染色体
s2.eG4 <- df[df$s2.eG4=="eG4"&df$chr1=="X",c(4,13)]
s2.eG4$type <- "s2.eG4"
colnames(s2.eG4) <- c("id","con","type")
s2.eG4$cell <- "S2"
kc.eG4 <- df[df$kc.eG4=="eG4"&df$chr1=="X",c(4,12)]
kc.eG4$type <- "Kc.eG4"
colnames(kc.eG4) <- c("id","con","type")
kc.eG4$cell <- "Kc167"
s2.noeG4 <- df[df$s2.eG4=="no eG4"&df$chr1=="X",c(4,13)]
s2.noeG4$type <- "S2noeG4"       
colnames(s2.noeG4) <- c("id","con","type")   
s2.noeG4$cell <- "S2"
kc.noeG4 <- df[df$kc.eG4=="no eG4"&df$chr1=="X",c(4,12)]
kc.noeG4$type <- "KcnoeG4"
colnames(kc.noeG4) <- c("id","con","type")
kc.noeG4$cell <- "Kc167"
df_combined <- rbind(s2.eG4, kc.eG4,s2.noeG4,kc.noeG4)
my_comparisons=list(c("KcnoeG4","Kc.eG4"),c("S2noeG4","s2.eG4"),c("Kc.eG4","s2.eG4"),c("KcnoeG4","S2noeG4"))  
df_combined$type <- factor(df_combined$type,levels = c("KcnoeG4","Kc.eG4","S2noeG4","s2.eG4"))
ggplot(data = df_combined,aes(x=type,y=con,fill=cell)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(-5000,-3000,-1800,-3800),
                     tip.length = 0,size=4) +
  scale_fill_manual(values = c("#4A6F8B","#CF4F39")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank()
  ) +
  ylab("Mean expression level") +
  coord_cartesian(ylim = c(-1000,15000)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title="DMSO(X chromosome)")+
  labs(fill = "Cell") +
  scale_x_discrete(labels = c("KcnoeG4" ="no eG4" , "S2noeG4" = "no eG4", "Kc.eG4" = "eG4", "s2.eG4" = "eG4")) 
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"DMSOKcS2.eG4noeG4.ExpressionLevel.ChrX.pdf"),
       device = "pdf",width = 4.5,height = 3.5)  

#可视化 常染色体和X染色体
s2.eG4 <- df[df$s2.eG4=="eG4",c(14,4,13)]
s2.eG4$type <- "s2.eG4"
colnames(s2.eG4) <- c("chr","id","con","type")
s2.eG4$cell <- "S2"
kc.eG4 <- df[df$kc.eG4=="eG4",c(14,4,12)]
kc.eG4$type <- "Kc.eG4"
colnames(kc.eG4) <- c("chr","id","con","type")
kc.eG4$cell <- "Kc167"
s2.noeG4 <- df[df$s2.eG4=="no eG4",c(14,4,13)]
s2.noeG4$type <- "S2noeG4"       
colnames(s2.noeG4) <- c("chr","id","con","type")
s2.noeG4$cell <- "S2"
kc.noeG4 <- df[df$kc.eG4=="no eG4",c(14,4,12)]
kc.noeG4$type <- "KcnoeG4"
colnames(kc.noeG4) <- c("chr","id","con","type")
kc.noeG4$cell <- "Kc167"
df_combined <- rbind(s2.eG4, kc.eG4,s2.noeG4,kc.noeG4)
my_comparisons=list(c("KcnoeG4","Kc.eG4"),c("S2noeG4","s2.eG4"),c("Kc.eG4","s2.eG4"),c("KcnoeG4","S2noeG4"))  
df_combined$type <- factor(df_combined$type,levels = c("KcnoeG4","Kc.eG4","S2noeG4","s2.eG4"))
ggplot(data = df_combined,aes(x=type,y=con,fill=cell)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(-42000,-39500,-38000,-40500),
                     tip.length = 0,size=4) +
  scale_fill_manual(values = c("#4A6F8B","#CF4F39")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank()
  ) +
  ylab("Mean expression level") +
  coord_cartesian(ylim = c(-1000,15000)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title="DMSO(All gene)")+
  labs(fill = "Cell") +
  scale_x_discrete(labels = c("KcnoeG4" ="no eG4" , "S2noeG4" = "no eG4", "Kc.eG4" = "eG4", "s2.eG4" = "eG4")) 
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"DMSOKcS2.eG4noeG4.ExpressionLevel.pdf"),
       device = "pdf",width = 7,height = 3.5)

#### 处理后 ####
rm(list = ls());gc();rm(list = ls())
Num = "010.7."
readcounts <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.counts.txt") %>% as.data.frame()
rownames(readcounts) <- readcounts$geneid
phen <- readcounts[,c(7:8,15:16)]
sample_name <- factor(colnames(phen))
metadata <- data.frame(sample_name)
metadata$cell <- as.factor(rep(c("kc","s2"),each=2))

dds <- DESeqDataSetFromMatrix(countData = phen,
                              colData = metadata,
                              design = ~ cell)
##2.dds标准化
dds <- DESeq(dds)
##3.获取标准化后的数据
normalized_counts <- counts(dds,normalized=T)
##4.有三个重复取平均
normalized_mean_counts = t(apply(normalized_counts, 1, function(a){tapply(a, metadata$cell, mean)}))

gene.kc.s2 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.2.gene.kc.s2.txt") %>% as.data.frame()
gene.kc.s2$kc.eG4 <- ifelse(gene.kc.s2$kc==1,"eG4","no eG4")
gene.kc.s2$s2.eG4 <- ifelse(gene.kc.s2$s2==2,"eG4","no eG4")
gene.kc.s2 <- subset(gene.kc.s2,chr=='2L'|chr=='2R'|chr=='3L'|chr=='3R'|chr=='X')

gene.kc.s2$kcphen <- normalized_mean_counts[match(gene.kc.s2$id,rownames(normalized_mean_counts)),1]
gene.kc.s2$s2phen <- normalized_mean_counts[match(gene.kc.s2$id,rownames(normalized_mean_counts)),2]

df <- na.omit(gene.kc.s2)
df$chr1 <- ifelse(df$chr=="X","X","A")
df$ratio <- df$kcphen/df$s2phen
df$ratio1 <- df$s2phen/df$kcphen
library(ggpubr)
my_comparisons = list(c("eG4","no eG4"))
ggplot(data = df,aes(x=kc.eG4,y=log2(ratio+1),fill=kc.eG4)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(4,4.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  # scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Control log2(Kc/S2)") +
  coord_cartesian(ylim = c(0,5)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title="Kc167 cell")

ggplot(data = df,aes(x=kc.eG4,y=log2(ratio+1),fill=kc.eG4)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  facet_grid(~chr1) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(4,4.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  # scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("PhenDC3 log2(Kc/S2+1)") +
  coord_cartesian(ylim = c(0,5)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title="Kc167 cell") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") 


df$ratio1 <- df$s2con/df$kccon
ggplot(data = df,aes(x=s2.eG4,y=log2(ratio1+1),fill=s2.eG4)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(4,4.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  # scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("PhenDC3 log2(S2/Kc+1)") +
  coord_cartesian(ylim = c(0,5)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title="S2 cell") 

ggplot(data = df,aes(x=s2.eG4,y=log2(ratio1+1),fill=s2.eG4)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  facet_grid(~chr1) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(4,4.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  # scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("PhenDC3 log2(S2/Kc+1)") +
  coord_cartesian(ylim = c(0,5)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title="S2 cell")+
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") 

# 比较两组独立样本的中位数是否存在差异。
wilcox.test(df[df$s2.eG4=="eG4"&df$chr1=="X",13],df[df$kc.eG4=="eG4"&df$chr1=="X",12]) #p-value = 0.0136
wilcox.test(df[df$s2.eG4=="no eG4"&df$chr1=="X",13],df[df$kc.eG4=="no eG4"&df$chr1=="X",12]) #0.3179
#可视化
s2.eG4 <- df[df$s2.eG4=="eG4"&df$chr1=="X",c(4,13)]
s2.eG4$type <- "s2.eG4"
colnames(s2.eG4) <- c("id","phen","type")
s2.eG4$cell <- "S2"
kc.eG4 <- df[df$kc.eG4=="eG4"&df$chr1=="X",c(4,12)]
kc.eG4$type <- "Kc.eG4"
colnames(kc.eG4) <- c("id","phen","type")
kc.eG4$cell <- "Kc167"
s2.noeG4 <- df[df$s2.eG4=="no eG4"&df$chr1=="X",c(4,13)]
s2.noeG4$type <- "S2noeG4"       
colnames(s2.noeG4) <- c("id","phen","type")   
s2.noeG4$cell <- "S2"
kc.noeG4 <- df[df$kc.eG4=="no eG4"&df$chr1=="X",c(4,12)]
kc.noeG4$type <- "KcnoeG4"
colnames(kc.noeG4) <- c("id","phen","type")
kc.noeG4$cell <- "Kc167"
df_combined <- rbind(s2.eG4, kc.eG4,s2.noeG4,kc.noeG4)
my_comparisons=list(c("KcnoeG4","Kc.eG4"),c("S2noeG4","s2.eG4"),c("Kc.eG4","s2.eG4"),c("KcnoeG4","S2noeG4"))  
df_combined$type <- factor(df_combined$type,levels = c("KcnoeG4","Kc.eG4","S2noeG4","s2.eG4"))
ggplot(data = df_combined,aes(x=type,y=phen,fill=cell)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(-7400,-4900,-3500,-6300),
                     tip.length = 0,size=4) +
 scale_fill_manual(values = c("#4A6F8B","#CF4F39")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank()
        ) +
  ylab("Mean expression level") +
  coord_cartesian(ylim = c(-1000,15000)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title="PhenDC3(X chromosome)")+
  labs(fill = "Cell") +
  scale_x_discrete(labels = c("KcnoeG4" ="no eG4" , "S2noeG4" = "no eG4", "Kc.eG4" = "eG4", "s2.eG4" = "eG4")) 
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"PhenKcS2.eG4noeG4.ExpressionLevel.ChrX.pdf"),
       device = "pdf",width = 4.5,height = 3.5)  

#可视化 常染色体和X染色体
s2.eG4 <- df[df$s2.eG4=="eG4",c(14,4,13)]
s2.eG4$type <- "s2.eG4"
colnames(s2.eG4) <- c("chr","id","con","type")
s2.eG4$cell <- "S2"
kc.eG4 <- df[df$kc.eG4=="eG4",c(14,4,12)]
kc.eG4$type <- "Kc.eG4"
colnames(kc.eG4) <- c("chr","id","con","type")
kc.eG4$cell <- "Kc167"
s2.noeG4 <- df[df$s2.eG4=="no eG4",c(14,4,13)]
s2.noeG4$type <- "S2noeG4"       
colnames(s2.noeG4) <- c("chr","id","con","type")
s2.noeG4$cell <- "S2"
kc.noeG4 <- df[df$kc.eG4=="no eG4",c(14,4,12)]
kc.noeG4$type <- "KcnoeG4"
colnames(kc.noeG4) <- c("chr","id","con","type")
kc.noeG4$cell <- "Kc167"
df_combined <- rbind(s2.eG4, kc.eG4,s2.noeG4,kc.noeG4)
my_comparisons=list(c("KcnoeG4","Kc.eG4"),c("S2noeG4","s2.eG4"),c("Kc.eG4","s2.eG4"),c("KcnoeG4","S2noeG4"))  
df_combined$type <- factor(df_combined$type,levels = c("KcnoeG4","Kc.eG4","S2noeG4","s2.eG4"))
ggplot(data = df_combined,aes(x=type,y=con,fill=cell)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(-44000,-41700,-40300,-42500),
                     tip.length = 0,size=4) +
  scale_fill_manual(values = c("#4A6F8B","#CF4F39")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank()
  ) +
  ylab("Mean expression level") +
  coord_cartesian(ylim = c(-1000,15000)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title="PhenDC3(All gene)")+
  labs(fill = "Cell") +
  scale_x_discrete(labels = c("KcnoeG4" ="no eG4" , "S2noeG4" = "no eG4", "Kc.eG4" = "eG4", "s2.eG4" = "eG4")) 
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"PhenKcS2.eG4noeG4.ExpressionLevel.pdf"),
       device = "pdf",width = 7,height = 3.5)

#### ####
Num = "010.6."
readcounts <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.counts.txt") %>% as.data.frame()
rownames(readcounts) <- readcounts$geneid
readcounts1 <- readcounts[,c(1:3,7:11,15:16)]
metadata <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.metadata.txt") %>% as.data.frame()
metadata1 <- metadata[c(1:3,7:11,15:16),]
metadata1$name <- as.factor(paste0(metadata1$treat,metadata1$cell))

dds <- DESeqDataSetFromMatrix(countData = readcounts1,
                              colData = metadata1,
                              design = ~ name)
##2.dds标准化
dds <- DESeq(dds)
##3.获取标准化后的数据
normalized_counts <- counts(dds,normalized=T)
##4.有三个重复取平均
normalized_mean_counts = t(apply(normalized_counts, 1, function(a){tapply(a, metadata1$name, mean)}))

kc.Phen <- na.omit(as.data.frame(results(dds, contrast = c("name", "Phenkc","conkc" ), cooksCutoff = FALSE))) #cooksCutoff = FALSE参数是指离群的基因不会用NA表示，之前没设置参数时，将会把离群值设为NA
kc.Phen$group <- ifelse(kc.Phen$pvalue<0.05&abs(kc.Phen$log2FoldChange)>=0.5,ifelse(kc.Phen$log2FoldChange>0.5,"Up","Down"),"No-sig")
table(kc.Phen$group)

s2.Phen <- na.omit(as.data.frame(results(dds, contrast = c("name", "Phens2","cons2" ), cooksCutoff = FALSE)))
s2.Phen$group <- ifelse(s2.Phen$pvalue<0.05&abs(s2.Phen$log2FoldChange)>=0.5,ifelse(s2.Phen$log2FoldChange>0.5,"Up","Down"),"No-sig")
table(s2.Phen$group)

gene.pqs <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.gene.pqs.bed") %>% as.data.frame()
colnames(gene.pqs) <- c("chr","start","end","id","strand","gene_symbol","gene.pqs")
gene.pqs$chr1 <- ifelse(gene.pqs$chr=="X","X","A")
gene.pqs$pqs.type <- ifelse(gene.pqs$gene.pqs==0,"other","pqs")
gene.kc.s2 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.2.gene.kc.s2.txt") %>% as.data.frame()
gene.pqs$kc <- gene.kc.s2[match(gene.pqs$id,gene.kc.s2$id),7] #匹配kc列
gene.pqs$kc.type <- ifelse(gene.pqs$pqs.type=="pqs",ifelse(gene.pqs$kc=="1","eG4","non-eG4"),"other") #增加kctype列
gene.pqs$s2 <- gene.kc.s2[match(gene.pqs$id,gene.kc.s2$id),8] #匹配s2列
gene.pqs$s2.type <- ifelse(gene.pqs$pqs.type=="pqs",ifelse(gene.pqs$s2=="2","eG4","non-eG4"),"other")

gene.pqs$conkc <- normalized_mean_counts[match(gene.pqs$id,rownames(normalized_mean_counts)),1]
gene.pqs$cons2 <- normalized_mean_counts[match(gene.pqs$id,rownames(normalized_mean_counts)),2]
gene.pqs$ratio <- gene.pqs$conkc/gene.pqs$cons2

dfkc <- subset(gene.pqs,chr=='2L'|chr=='2R'|chr=='3L'|chr=='3R'|chr=='X')
dfkc$group <- kc.Phen[match(dfkc$id,rownames(kc.Phen)),7]
dfkc$type <- paste0(dfkc$kc.type,"_",dfkc$group)
dfkc <- na.omit(dfkc)
non_finite_rows <- which(!is.finite(dfkc$ratio))

# 如果需要，可以将非有限值替换为 NA
dfkc$ratio[non_finite_rows] <- NA

# 删除包含缺失值的行
dfkc_clean <- dfkc[complete.cases(dfkc), ]
dfkc_clean$type <- factor(dfkc_clean$type,levels = c("eG4_Up","eG4_Down","eG4_No-sig","non-eG4_Up","non-eG4_Down",
                                         "non-eG4_No-sig","other_Up","other_Down","other_No-sig"))
my_comparisons = list(c("eG4_Up","eG4_Down"),c("eG4_Up","eG4_No-sig"),c("eG4_Down","eG4_No-sig"),
                      c("non-eG4_Up","non-eG4_Down"),c("non-eG4_Up","non-eG4_No-sig"),c("non-eG4_Down","non-eG4_No-sig"),
                      c("other_Up","other_Down"),c("other_Up","other_No-sig"),c("other_Down","other_No-sig"))

gene_counts <- as.data.frame(table(dfkc_clean$type))
colnames(gene_counts) <- c("type", "gene_count")

ggplot(data = dfkc_clean, aes(x = type, y = log2(ratio + 1), fill = type)) +
  geom_boxplot(notch = FALSE, outlier.colour = "white") +
  cowplot::theme_half_open() +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = rep(c(4, 4.3, 4.6), 3),
                     aes(label = paste0("p = ", ..p.format..)), tip.length = 0, size = 4) +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Control log2(Kc/S2)") +
  coord_cartesian(ylim = c(0, 6)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Kc167 cell") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_text(data = gene_counts, aes(label = gene_count, y = 3), vjust = 0, size = 4) +
  geom_text(data = gene_counts, aes(label = type, y = 3.5), vjust = 0, size = 4, color = "red")

ggplot(data = dfkc_clean, aes(x = type, y = log2(ratio + 1), fill = type)) +
  geom_boxplot(notch = FALSE, outlier.colour = "white") +
  facet_grid(~chr1) +
  cowplot::theme_half_open() +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = rep(c(4, 4.3, 4.6), 3),
                     aes(label = paste0("p = ", ..p.format..)), tip.length = 0, size = 4) +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Control log2(Kc/S2)") +
  coord_cartesian(ylim = c(0, 6)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Kc167 cell") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_text(data = gene_counts, aes(label = gene_count, y = 3), vjust = 0, size = 4) +
  geom_text(data = gene_counts, aes(label = type, y = 3.5), vjust = 0, size = 4, color = "red")

my_comparisons = list(c("eG4_Up","non-eG4_Up"),c("eG4_Up","other_Up"),c("non-eG4_Up","other_Up"),
                      c("eG4_Down","non-eG4_Down"),c("eG4_Down","other_Down"),c("non-eG4_Down","other_Down"),
                      c("eG4_No-sig","non-eG4_No-sig"),c("eG4_No-sig","other_No-sig"),c("non-eG4_No-sig","other_No-sig"))
ggplot(data = dfkc_clean, aes(x = type, y = log2(ratio + 1), fill = type)) +
  geom_boxplot(notch = FALSE, outlier.colour = "white") +
  cowplot::theme_half_open() +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(3.7,4, 4.3, 4.6,4.9,5.2,5.5,5.8,6.1,6.4),
                     aes(label = paste0("p = ", ..p.format..)), tip.length = 0, size = 4) +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Control log2(Kc/S2)") +
  coord_cartesian(ylim = c(0, 8)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Kc167 cell") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_text(data = gene_counts, aes(label = gene_count, y = 3), vjust = 0, size = 4) +
  geom_text(data = gene_counts, aes(label = type, y = 3.5), vjust = 0, size = 4, color = "red")

gene_counts <- as.data.frame(table(dfkc_clean$type1))
colnames(gene_counts) <- c("type1", "gene_count")
ggplot(data = dfkc_clean, aes(x = type, y = log2(ratio + 1), fill = type)) +
  geom_boxplot(notch = FALSE, outlier.colour = "white") +
  facet_grid(~chr1)+
  cowplot::theme_half_open() +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(3.7,4, 4.3, 4.6,4.9,5.2,5.5,5.8,6.1,6.4),
                     aes(label = paste0("p = ", ..p.format..)), tip.length = 0, size = 4) +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank()) +
  ylab("Control log2(Kc/S2)") +
  coord_cartesian(ylim = c(0, 8)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Kc167 cell") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") 




gene.pqs <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.gene.pqs.bed") %>% as.data.frame()
colnames(gene.pqs) <- c("chr","start","end","id","strand","gene_symbol","gene.pqs")
gene.pqs$chr1 <- ifelse(gene.pqs$chr=="X","X","A")
gene.pqs$pqs.type <- ifelse(gene.pqs$gene.pqs==0,"other","pqs")
gene.kc.s2 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.2.gene.kc.s2.txt") %>% as.data.frame()
gene.pqs$kc <- gene.kc.s2[match(gene.pqs$id,gene.kc.s2$id),7] #匹配kc列
gene.pqs$kc.type <- ifelse(gene.pqs$pqs.type=="pqs",ifelse(gene.pqs$kc=="1","eG4","non-eG4"),"other") #增加kctype列
gene.pqs$s2 <- gene.kc.s2[match(gene.pqs$id,gene.kc.s2$id),8] #匹配s2列
gene.pqs$s2.type <- ifelse(gene.pqs$pqs.type=="pqs",ifelse(gene.pqs$s2=="2","eG4","non-eG4"),"other")

gene.pqs$conkc <- normalized_mean_counts[match(gene.pqs$id,rownames(normalized_mean_counts)),1]
gene.pqs$cons2 <- normalized_mean_counts[match(gene.pqs$id,rownames(normalized_mean_counts)),2]
gene.pqs$ratio <- gene.pqs$cons2/gene.pqs$conkc

dfs2 <- subset(gene.pqs,chr=='2L'|chr=='2R'|chr=='3L'|chr=='3R'|chr=='X')
dfs2$group <- s2.Phen[match(dfs2$id,rownames(s2.Phen)),7]
dfs2$type <- paste0(dfs2$s2.type,"_",dfs2$group)
dfs2 <- na.omit(dfs2)

dfs2$type <- factor(dfs2$type,levels = c("eG4_Up","eG4_Down","eG4_No-sig","non-eG4_Up","non-eG4_Down",
                                                     "non-eG4_No-sig","other_Up","other_Down","other_No-sig"))
my_comparisons = list(c("eG4_Up","eG4_Down"),c("eG4_Up","eG4_No-sig"),c("eG4_Down","eG4_No-sig"),
                      c("non-eG4_Up","non-eG4_Down"),c("non-eG4_Up","non-eG4_No-sig"),c("non-eG4_Down","non-eG4_No-sig"),
                      c("other_Up","other_Down"),c("other_Up","other_No-sig"),c("other_Down","other_No-sig"))

gene_counts <- as.data.frame(table(dfs2$type))
colnames(gene_counts) <- c("type", "gene_count")

ggplot(data = dfs2, aes(x = type, y = log2(ratio + 1), fill = type)) +
  geom_boxplot(notch = FALSE, outlier.colour = "white") +
  cowplot::theme_half_open() +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = rep(c(4, 4.3, 4.6), 3),
                     aes(label = paste0("p = ", ..p.format..)), tip.length = 0, size = 4) +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Control log2(S2/Kc)") +
  coord_cartesian(ylim = c(0, 6)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "S2 cell") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_text(data = gene_counts, aes(label = gene_count, y = 3), vjust = 0, size = 4) +
  geom_text(data = gene_counts, aes(label = type, y = 3.5), vjust = 0, size = 4, color = "red")

ggplot(data = dfs2, aes(x = type, y = log2(ratio + 1), fill = type)) +
  geom_boxplot(notch = FALSE, outlier.colour = "white") +
  facet_grid(~chr1) +
  cowplot::theme_half_open() +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = rep(c(4, 4.3, 4.6), 3),
                     aes(label = paste0("p = ", ..p.format..)), tip.length = 0, size = 4) +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Control log2(S2/Kc)") +
  coord_cartesian(ylim = c(0, 6)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "S2 cell") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") 


my_comparisons = list(c("eG4_Up","non-eG4_Up"),c("eG4_Up","other_Up"),c("non-eG4_Up","other_Up"),
                      c("eG4_Down","non-eG4_Down"),c("eG4_Down","other_Down"),c("non-eG4_Down","other_Down"),
                      c("eG4_No-sig","non-eG4_No-sig"),c("eG4_No-sig","other_No-sig"),c("non-eG4_No-sig","other_No-sig"))
ggplot(data = dfs2, aes(x = type, y = log2(ratio + 1), fill = type)) +
  geom_boxplot(notch = FALSE, outlier.colour = "white") +
  cowplot::theme_half_open() +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(3.7,4, 4.3, 4.6,4.9,5.2,5.5,5.8,6.1,6.4),
                     aes(label = paste0("p = ", ..p.format..)), tip.length = 0, size = 4) +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Control log2(S2/Kc)") +
  coord_cartesian(ylim = c(0, 8)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "S2 cell") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_text(data = gene_counts, aes(label = gene_count, y = 3), vjust = 0, size = 4) +
  geom_text(data = gene_counts, aes(label = type, y = 3.5), vjust = 0, size = 4, color = "red")

ggplot(data = dfs2, aes(x = type, y = log2(ratio + 1), fill = type)) +
  geom_boxplot(notch = FALSE, outlier.colour = "white") +
  facet_grid(~chr1) +
  cowplot::theme_half_open() +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(3.7,4, 4.3, 4.6,4.9,5.2,5.5,5.8,6.1,6.4),
                     aes(label = paste0("p = ", ..p.format..)), tip.length = 0, size = 4) +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Control log2(S2/Kc)") +
  coord_cartesian(ylim = c(0, 8)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "S2 cell") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") 

dfs2$type1 <- paste(dfs2$chr1,dfs2$type)
gene_counts <- as.data.frame(table(dfs2$type1))
colnames(gene_counts) <- c("type1", "gene_count")
