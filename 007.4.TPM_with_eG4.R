## bedtools intersect -a "/home/yuss/flyG4/data/ref/gene.sort.bed" -b /home/yuss/flyG4/result/PQS/001.2.dmel.pqs.bed -c > /home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.gene.pqs.bed

#*基因上含有不同eG4的TPM-------------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())#清空
Num = "007.4."
#### 计算TPM ####
#之前计算过基因长度
gene.length <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.2.dmel.genelength.txt") %>% as.data.frame()
counts <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.1.counts.txt") %>% as.data.frame()
tpm.calculate = function(exprset,len){
  readperlength = t(do.call(rbind, lapply(1:ncol(exprset), function(i){
    exprset[,i]/len})))
  totalcounts <- colSums(readperlength)
  tpm = t(apply(readperlength, 1, function(x) 10^6 * x/totalcounts)) %>% as.data.frame()
  colnames(tpm) = colnames(exprset)
  row.names(tpm) = row.names(exprset)
  return(tpm)
}
counts$length = gene.length[match(counts$geneid,gene.length$gene.name),2]
counts <- column_to_rownames(counts,var= "geneid")
tpm = tpm.calculate(counts[,-5],counts$length) #ncol() 函数返回矩阵的列数
tpm$aver_Kc <- rowMeans(tpm[,c("Kc167-2","Kc167-4")])
tpm$aver_S2 <- rowMeans(tpm[,c("S2-DRSC-14","S2-DRSC-16")])
tpm$gene_id <- rownames(tpm)
write.table(tpm,file = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/",Num,"Tpm.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)
#### 分成三类other,eG4,non-eG4 ####
gene.pqs <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.gene.pqs.bed") %>% as.data.frame()
colnames(gene.pqs) <- c("chr","start","end","id","strand","gene_symbol","gene.pqs")
gene.pqs$pqs.type <- ifelse(gene.pqs$gene.pqs==0,"other","pqs")
gene.kc.s2 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.2.gene.kc.s2.txt") %>% as.data.frame()
gene.pqs$kc <- gene.kc.s2[match(gene.pqs$id,gene.kc.s2$id),7] #匹配kc列
gene.pqs$kc.type <- ifelse(gene.pqs$pqs.type=="pqs",ifelse(gene.pqs$kc=="1","eG4","non-eG4"),"other") #增加kctype列
gene.pqs$s2 <- gene.kc.s2[match(gene.pqs$id,gene.kc.s2$id),8] #匹配s2列
gene.pqs$s2.type <- ifelse(gene.pqs$pqs.type=="pqs",ifelse(gene.pqs$s2=="2","eG4","non-eG4"),"other")
#把tpm匹配到gene.pqs，tpm17494,gene17754,匹配不到的表达量设为0
tpm <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.Tpm.txt") %>% as.data.frame()
gene.pqs$kc.tpm <- tpm[match(gene.pqs$id,tpm$gene_id),5]
# gene.pqs[is.na(gene.pqs$kc.tpm), 13]=0
gene.pqs$s2.tpm <- tpm[match(gene.pqs$id,tpm$gene_id),6]
# gene.pqs[is.na(gene.pqs$s2.tpm), 14]=0
#另一种
gene.pqs <- subset(gene.pqs,s2.tpm!='NA')
#画图
gene.pqs$kc.type <- factor(gene.pqs$kc.type,levels = c("non-eG4","eG4","other"))
gene.pqs$s2.type <- factor(gene.pqs$s2.type,levels = c("non-eG4","eG4","other"))
my_comparisons = list(c("non-eG4","eG4"),c("eG4","other"))
ggplot(data = gene.pqs,aes(x=kc.type,y=log2(kc.tpm+1),fill=kc.type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(12.5,14),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="All genes (Kc167 cells)") +
  coord_cartesian(ylim = c(0,15)) 
#ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"eG4_TPM_in_Kc.pdf"),
       device = "pdf",width = 3.5,height = 3)
ggplot(data = gene.pqs,aes(x=s2.type,y=log2(s2.tpm+1),fill=s2.type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(12.5,13.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#4d4d4d","#D7A49D","#92c5de")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="All genes (S2 cells)") +
  coord_cartesian(ylim = c(0,15)) 
#ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"eG4_TPM_in_S2.pdf"),
         device = "pdf",width = 3.5,height = 3)

aggregate(s2.tpm ~ s2.type, data = gene.pqs, FUN = median)

#### 统计含有eG4基因表达的比例 ####
gene.pqs$kc.tpm.type <- ifelse(gene.pqs$kc.tpm>=1,"express","nonexpress")
table(gene.pqs$kc.tpm.type,gene.pqs$kc.type)
df <- table(gene.pqs$kc.tpm.type,gene.pqs$kc.type) %>% as.data.frame()
#长数据转为宽数据
df <- spread(df,Var1,Freq)
df$sum <- rowSums(df[,2:3])
df$express_ratio <- df$express/df$sum

gene.pqs$s2.tpm.type <- ifelse(gene.pqs$s2.tpm>1,"express","nonexpress")
df <- table(gene.pqs$s2.tpm.type,gene.pqs$s2.type) %>% as.data.frame()
#长数据转为宽数据
df <- spread(df,Var1,Freq)
df$sum <- rowSums(df[,2:3])
df$express_ratio <- df$express/df$sum

#### 分成两类no eG4,eG4 ####
rm(list = ls());gc();rm(list = ls())
gene.kc.s2 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.2.gene.kc.s2.txt") %>% as.data.frame()
gene.kc.s2$kc.eG4 <- ifelse(gene.kc.s2$kc==1,"eG4","no eG4")
gene.kc.s2$s2.eG4 <- ifelse(gene.kc.s2$s2==2,"eG4","no eG4")
#把tpm匹配到gene.pqs，tpm17494,gene17754,匹配不到的表达量设为0
tpm <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.Tpm.txt") %>% as.data.frame()
gene.kc.s2$kc.tpm <- tpm[match(gene.kc.s2$id,tpm$gene_id),5]
# gene.kc.s2[is.na(gene.kc.s2$kc.tpm), 12]=0 ##设置后在kc和s2细胞中会促进基因表达
gene.kc.s2$s2.tpm <- tpm[match(gene.kc.s2$id,tpm$gene_id),6]
# gene.kc.s2[is.na(gene.kc.s2$s2.tpm), 13]=0
#另一种
gene.kc.s2 <- subset(gene.kc.s2,s2.tpm!='NA')
#画图
gene.kc.s2$kc.eG4 <- factor(gene.kc.s2$kc.eG4,levels = c("no eG4","eG4"))
gene.kc.s2$s2.eG4 <- factor(gene.kc.s2$s2.eG4,levels = c("no eG4","eG4"))
my_comparisons = list(c("no eG4","eG4"))
ggplot(data = gene.kc.s2,aes(x=kc.eG4,y=log2(kc.tpm+1),fill=kc.eG4)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     method="t.test",
                     paired = FALSE,  # 设置为 FALSE，表示独立样本 t 检验
                     label.y = c(12.5,14),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="All genes (Kc167 cells)") +
  coord_cartesian(ylim = c(0,15)) 
#ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"eG4_TPM_in_Kc.pdf"),
device = "pdf",width = 3.5,height = 3)
library(ggpubr)
ggplot(data = gene.kc.s2,aes(x=s2.eG4,y=log2(s2.tpm+1),fill=s2.eG4)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     method="t.test",
                     paired = FALSE,  # 设置为 FALSE，表示独立样本 t 检验
                     label.y = c(12.5,13.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#4d4d4d","#D7A49D","#92c5de")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="All genes (S2 cells)") +
  coord_cartesian(ylim = c(0,15)) 
#ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"eG4_TPM_in_S2.pdf"),
device = "pdf",width = 3.5,height = 3)



aggregate(s2.tpm ~ s2.eG4, data = gene.kc.s2, FUN = median)
# > median
# s2.eG4   s2.tpm
# 1 no eG4 3.773099
# 2    eG4 3.776719
#比较平均数
s2.noeG4 <- gene.kc.s2 %>% filter(s2.eG4=="no eG4")
s2.eG4 <- gene.kc.s2 %>% filter(s2.eG4=="eG4")
t.test(s2.noeG4$s2.tpm,s2.eG4$s2.tpm)
# t = 2.9151, df = 7251.3, p-value = 0.003567
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   8.714199 44.497034
# sample estimates:
#   mean of x mean of y 
# 80.47995  53.87434


#### 基因组区域两类G4在不同染色体上(5个)的表达水平 ####
subset_gene.kc.s2 <- subset(gene.kc.s2, chr %in% c("2L", "2R", "3L","3R","X"))
my_comparisons = list(c("no eG4","eG4"))

ggplot(data = subset_gene.kc.s2,aes(x=kc.eG4,y=log2(kc.tpm+1),fill=kc.eG4)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  facet_grid(~ chr) +
  stat_compare_means(comparisons = my_comparisons,
                     method="t.test",
                     paired = FALSE,  # 设置为 FALSE，表示独立样本 t 检验
                     label.y = rep(13,3),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#4d4d4d","#FDDBC7")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="All genes (Kc167 cells)") +
  coord_cartesian(ylim = c(0,15)) 
table(subset_gene.kc.s2$chr,subset_gene.kc.s2$kc.eG4)
aggregate(kc.tpm ~ kc.eG4 + chr , data = subset_gene.kc.s2, FUN = median)
# kc.eG4 chr   kc.tpm
# 1  no eG4  2L 4.103364
# 2     eG4  2L 1.346166
# 3  no eG4  2R 6.449765
# 4     eG4  2R 6.657282
# 5  no eG4  3L 5.392459
# 6     eG4  3L 3.479306
# 7  no eG4  3R 5.596446
# 8     eG4  3R 4.514403
# 9  no eG4   X 9.122236
# 10    eG4   X 5.639902


ggplot(data = subset_gene.kc.s2,aes(x=s2.eG4,y=log2(s2.tpm+1),fill=s2.eG4)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  facet_grid(~ chr) +
  stat_compare_means(comparisons = my_comparisons,
                     method="t.test",
                     paired = FALSE,  # 设置为 FALSE，表示独立样本 t 检验
                     label.y = rep(13,3),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#4d4d4d","#D7A49D")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="All genes (S2 cells)") +
  coord_cartesian(ylim = c(0,15))
aggregate(s2.tpm ~ s2.eG4 + chr , data = subset_gene.kc.s2, FUN = mean)

a <- subset_gene.kc.s2[subset_gene.kc.s2$chr=="X",]
a1 <- spread(a,s2.eG4,s2.tpm) ##长数据转换为宽数据，Var2为需要分解的变量，Freq为分解后的列的取值
t.test(a1$`no eG4`, a1$eG4)

#### X染色体上含有G4的基因的表达比例 ####
X.gene.kc.s2 <- gene.kc.s2 %>% filter(chr=="X")
X.gene.kc.s2$kc.tpm.type <- ifelse(X.gene.kc.s2$kc.tpm>=1,"express","nonexpress")
table(X.gene.kc.s2$kc.tpm.type,X.gene.kc.s2$kc.eG4)
kc.df <- table(X.gene.kc.s2$kc.tpm.type,X.gene.kc.s2$kc.eG4) %>% as.data.frame()
#长数据转为宽数据
kc.df <- spread(kc.df,Var1,Freq)
kc.df$sum <- rowSums(kc.df[,2:3])
kc.df$express_ratio <- kc.df$express/kc.df$sum
kc.s2.express <- kc.df[,c(1,5)]
X.gene.kc.s2$s2.tpm.type <- ifelse(X.gene.kc.s2$s2.tpm>=1,"express","nonexpress")
table(X.gene.kc.s2$s2.tpm.type,X.gene.kc.s2$s2.eG4)
s2.df <- table(X.gene.kc.s2$s2.tpm.type,X.gene.kc.s2$s2.eG4) %>% as.data.frame()
#长数据转为宽数据
s2.df <- spread(s2.df,Var1,Freq)
s2.df$sum <- rowSums(s2.df[,2:3])
s2.df$express_ratio <- s2.df$express/s2.df$sum
kc.s2.express$s2 <- s2.df[,5]
colnames(kc.s2.express)[2] <- "kc"
kc.s2.express.long <- pivot_longer(kc.s2.express, cols = c(kc, s2), names_to = "Measure", values_to = "Value")
kc.s2.express.long <- kc.s2.express.long %>% mutate(type = paste(Var2, Measure,sep = "_"))
kc.s2.express.long$type <- factor(kc.s2.express.long$type,levels = c("no eG4_kc","eG4_kc","no eG4_s2","eG4_s2"))
kc.s2.express.long$Value <- round(kc.s2.express.long$Value,2)*100
kc.s2.express.long$percent <- paste(kc.s2.express.long$Value, "%",sep = "")

ggplot(kc.s2.express.long,aes(x=type,y=Value,fill=type))+
  geom_bar(stat = "identity",position = "dodge",width = 0.8)+
  theme_bw()+ 
  ylab("% Expressed genes")+
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=10,colour = "black"),axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,colour = "black"),axis.title.y = element_text(size=12),
        legend.position = "none")+
  scale_y_continuous(expand = c(0,0),limits = c(0,70)) +
  geom_vline(xintercept = 2.5,color="#595959",size=0.4) +
  geom_text(aes(label=percent,vjust=-0.5),color="black",size=4) + 
  labs(title="Genes on the X chromosome (Kc167 and S2 cells)") +
  scale_fill_manual(values = c("no eG4_kc"="#8EA1CC","eG4_kc"="#FC8E62","no eG4_s2"="#8EA1CC","eG4_s2"="#66C3AA")) +
  scale_x_discrete(labels = c("no eG4","eG4","no eG4","eG4"))


#### G4密度 ####
subset_gene.kc.s2$length <- subset_gene.kc.s2$end-subset_gene.kc.s2$start
table(subset_gene.kc.s2$length,subset_gene.kc.s2$chr)
sum(subset_gene.kc.s2$chr=="2L",subset_gene.kc.s2$length)
sum(subset_gene.kc.s2$chr=="2L",subset_gene.kc.s2$kc)
for (i in unique(subset_gene.kc.s2$chr)){
  print(i)
  # print(sum(subset_gene.kc.s2$chr==i,subset_gene.kc.s2$kc))
  # print(sum(subset_gene.kc.s2$chr==i,subset_gene.kc.s2$length))
  print(sum(subset_gene.kc.s2$chr==i,subset_gene.kc.s2$kc)/sum(subset_gene.kc.s2$chr==i,subset_gene.kc.s2$length))}

for (i in unique(subset_gene.kc.s2$chr)){
  print(i)
  # print(sum(subset_gene.kc.s2$chr==i,subset_gene.kc.s2$kc))
  # print(sum(subset_gene.kc.s2$chr==i,subset_gene.kc.s2$length))
  print(sum(subset_gene.kc.s2$chr==i,subset_gene.kc.s2$s2)/sum(subset_gene.kc.s2$chr==i,subset_gene.kc.s2$length))}

#### 基因组区域两类G4在常染色体（A）、性染色体上的表达水平 ####
subset_gene.kc.s2 <- subset_gene.kc.s2 %>%
  mutate(chr_type = ifelse(chr %in% c("2L", "2R", "3L", "3R"), "A", "X"))
my_comparisons <- list(c("no eG4","eG4"))
ggplot(data = subset_gene.kc.s2,aes(x=kc.eG4,y=log2(kc.tpm+1),fill=kc.eG4)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  facet_grid(~ chr_type) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(13,13),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#4d4d4d","#FDDBC7")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="All genes (Kc167 cells)") +
  coord_cartesian(ylim = c(0,15)) 

ggplot(data = subset_gene.kc.s2,aes(x=s2.eG4,y=log2(s2.tpm+1),fill=s2.eG4)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  facet_grid(~ chr_type) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(13,13),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#4d4d4d","#D7A49D")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="All genes (S2 cells)") +
  coord_cartesian(ylim = c(0,15))
sum(subset_gene.kc.s2$length,subset_gene.kc.s2$chr)

#### 在A和X染色体上含有G4基因的比例 ####
# table(subset_gene.kc.s2$chr,subset_gene.kc.s2$kc.eG4, subset_gene.kc.s2$kc.tpm.type)
# kc.chr.ratio <- table(subset_gene.kc.s2$chr,subset_gene.kc.s2$kc.eG4, subset_gene.kc.s2$kc.tpm.type) %>% as.data.frame()
df <- table(subset_gene.kc.s2$kc.eG4, subset_gene.kc.s2$chr_type) %>% as.data.frame()
df <- spread(df,Var1,Freq)
df$sum <- rowSums(df[,2:3])
df$containing_ratio <- df$eG4/df$sum
kc.s2.containing <- df[,c(1,5)]

df <- table(subset_gene.kc.s2$s2.eG4, subset_gene.kc.s2$chr_type) %>% as.data.frame()
df <- spread(df,Var1,Freq)
df$sum <- rowSums(df[,2:3])
df$containing_ratio <- df$eG4/df$sum
kc.s2.containing$s2 <- df[,5]
colnames(kc.s2.containing)[2] <- "kc"
kc.s2.containing.long <- pivot_longer(kc.s2.containing, cols = c(kc, s2), names_to = "Measure", values_to = "Value")
kc.s2.containing.long <- kc.s2.containing.long %>% mutate(type = paste(Var2, Measure,sep = "_"))
kc.s2.containing.long$type <- factor(kc.s2.containing.long$type,levels = c("A_kc","X_kc","A_s2","X_s2"))
kc.s2.containing.long$Value <- round(kc.s2.containing.long$Value,2)*100
kc.s2.containing.long$percent <- paste(kc.s2.containing.long$Value, "%",sep = "")
ggplot(kc.s2.containing.long,aes(x=type,y=Value,fill=type))+
  geom_bar(stat = "identity",position = "dodge",width = 0.8)+
  theme_bw()+ 
  ylab("% eG4 containing genes")+
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=10,colour = "black"),axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,colour = "black"),axis.title.y = element_text(size=12),
        legend.position = "none")+
  scale_y_continuous(expand = c(0,0),limits = c(0,30)) +
  geom_vline(xintercept = 2.5,color="#595959",size=0.4) +
  geom_text(aes(label=percent,vjust=-0.5),color="black",size=4) + 
  labs(title="All genes (Kc167 and S2 cells)") +
  scale_fill_manual(values = c("#8EA1CC","#FC8E62","#8EA1CC","#66C3AA"))# +
#  scale_x_discrete(labels = c("A-Kc","eG4","no eG4","eG4"))
#ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"Promoter.noeG4.eG4.ExpressRatio_in_KcS2.pdf"),
       device = "pdf",width = 4,height = 4.8)


#### 统计含有eG4基因表达的比例 ####
gene.kc.s2$kc.tpm.type <- ifelse(gene.kc.s2$kc.tpm>=1,"express","nonexpress")
table(gene.kc.s2$kc.tpm.type,gene.kc.s2$kc.eG4)
kc.df <- table(gene.kc.s2$kc.tpm.type,gene.kc.s2$kc.eG4) %>% as.data.frame()
#长数据转为宽数据
kc.df <- spread(kc.df,Var1,Freq)
kc.df$sum <- rowSums(kc.df[,2:3])
kc.df$express_ratio <- kc.df$express/kc.df$sum
kc.s2.express <- kc.df[,c(1,5)]
gene.kc.s2$s2.tpm.type <- ifelse(gene.kc.s2$s2.tpm>=1,"express","nonexpress")
table(gene.kc.s2$s2.tpm.type,gene.kc.s2$s2.eG4)
s2.df <- table(gene.kc.s2$s2.tpm.type,gene.kc.s2$s2.eG4) %>% as.data.frame()
#长数据转为宽数据
s2.df <- spread(s2.df,Var1,Freq)
s2.df$sum <- rowSums(s2.df[,2:3])
s2.df$express_ratio <- s2.df$express/s2.df$sum
kc.s2.express$s2 <- s2.df[,5]
colnames(kc.s2.express)[2] <- "kc"
kc.s2.express.long <- pivot_longer(kc.s2.express, cols = c(kc, s2), names_to = "Measure", values_to = "Value")
kc.s2.express.long <- kc.s2.express.long %>% mutate(type = paste(Var2, Measure,sep = "_"))
kc.s2.express.long$type <- factor(kc.s2.express.long$type,levels = c("no eG4_kc","eG4_kc","no eG4_s2","eG4_s2"))
kc.s2.express.long$Value <- round(kc.s2.express.long$Value,2)*100
kc.s2.express.long$percent <- paste(kc.s2.express.long$Value, "%",sep = "")
ggplot(kc.s2.express.long,aes(x=type,y=Value,fill=type))+
  geom_bar(stat = "identity",position = "dodge",width = 0.8)+
  theme_bw()+ 
  ylab("% Expressed genes")+
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=10,colour = "black"),axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,colour = "black"),axis.title.y = element_text(size=12),
        legend.position = "none")+
  scale_y_continuous(expand = c(0,0),limits = c(0,70)) +
  geom_vline(xintercept = 2.5,color="#595959",size=0.4) +
  geom_text(aes(label=percent,vjust=-0.5),color="black",size=4) + 
  labs(title="All genes (Kc167 and S2 cells)") +
  scale_fill_manual(values = c("no eG4_kc"="#8EA1CC","eG4_kc"="#FC8E62","no eG4_s2"="#8EA1CC","eG4_s2"="#66C3AA")) +
  scale_x_discrete(labels = c("no eG4","eG4","no eG4","eG4"))
#ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"Promoter.noeG4.eG4.ExpressRatio_in_KcS2.pdf"),
       device = "pdf",width = 4,height = 4.8)
## sed '1d' /home/yuss/flyG4/data/ref/promoter.bed | bedtools intersect -a - -b /home/yuss/flyG4/result/PQS/001.2.dmel.pqs.bed -c > /home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter.pqs.bed


#### 排除了4 Y染色体 ####
subset_gene.kc.s2$kc.tpm.type <- ifelse(subset_gene.kc.s2$kc.tpm>=1,"express","nonexpress")
table(subset_gene.kc.s2$kc.tpm.type,subset_gene.kc.s2$kc.eG4)
kc.df <- table(subset_gene.kc.s2$kc.tpm.type,subset_gene.kc.s2$kc.eG4) %>% as.data.frame()
#长数据转为宽数据
kc.df <- spread(kc.df,Var1,Freq)
kc.df$sum <- rowSums(kc.df[,2:3])
kc.df$express_ratio <- kc.df$express/kc.df$sum
kc.s2.express <- kc.df[,c(1,5)]
subset_gene.kc.s2$s2.tpm.type <- ifelse(subset_gene.kc.s2$s2.tpm>=1,"express","nonexpress")
table(subset_gene.kc.s2$s2.tpm.type,subset_gene.kc.s2$s2.eG4)
s2.df <- table(subset_gene.kc.s2$s2.tpm.type,subset_gene.kc.s2$s2.eG4) %>% as.data.frame()
#长数据转为宽数据
s2.df <- spread(s2.df,Var1,Freq)
s2.df$sum <- rowSums(s2.df[,2:3])
s2.df$express_ratio <- s2.df$express/s2.df$sum
kc.s2.express$s2 <- s2.df[,5]
colnames(kc.s2.express)[2] <- "kc"
kc.s2.express.long <- pivot_longer(kc.s2.express, cols = c(kc, s2), names_to = "Measure", values_to = "Value")
kc.s2.express.long <- kc.s2.express.long %>% mutate(type = paste(Var2, Measure,sep = "_"))
kc.s2.express.long$type <- factor(kc.s2.express.long$type,levels = c("no eG4_kc","eG4_kc","no eG4_s2","eG4_s2"))
kc.s2.express.long$Value <- round(kc.s2.express.long$Value,2)*100
kc.s2.express.long$percent <- paste(kc.s2.express.long$Value, "%",sep = "")

ggplot(kc.s2.express.long,aes(x=type,y=Value,fill=type))+
  geom_bar(stat = "identity",position = "dodge",width = 0.8)+
  theme_bw()+ 
  ylab("% Expressed genes")+
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=10,colour = "black"),axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,colour = "black"),axis.title.y = element_text(size=12),
        legend.position = "none")+
  scale_y_continuous(expand = c(0,0),limits = c(0,70)) +
  geom_vline(xintercept = 2.5,color="#595959",size=0.4) +
  geom_text(aes(label=percent,vjust=-0.5),color="black",size=4) + 
  labs(title="All genes (Kc167 and S2 cells)") +
  scale_fill_manual(values = c("no eG4_kc"="#8EA1CC","eG4_kc"="#FC8E62","no eG4_s2"="#8EA1CC","eG4_s2"="#66C3AA")) +
  scale_x_discrete(labels = c("no eG4","eG4","no eG4","eG4"))


#*promoter区：基因上下游1500bp 基因含有eG4和不含有eG4的TPM-----------------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())#清空
Num = "007.4."
#### 导入TPM ####
tpm <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.Tpm.txt") %>% as.data.frame()
#### 分成三类other,eG4,non-eG4 ####
promoter.pqs <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter.pqs.bed") %>% as.data.frame()
colnames(promoter.pqs) <- c("chr","start","end","id","strand","promoter.pqs")
promoter.pqs$pqs.type <- ifelse(promoter.pqs$promoter.pqs==0,"other","pqs")
promoter.kc.s2 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.2.promoter.kc.s2.txt") %>% as.data.frame()
promoter.pqs$kc <- promoter.kc.s2[match(promoter.pqs$id,promoter.kc.s2$gene_id),6] #匹配kc列
promoter.pqs$kc.type <- ifelse(promoter.pqs$pqs.type=="pqs",ifelse(promoter.pqs$kc=="1","eG4","non-eG4"),"other") #增加kctype列
promoter.pqs$s2 <- promoter.kc.s2[match(promoter.pqs$id,promoter.kc.s2$gene_id),7] #匹配s2列
promoter.pqs$s2.type <- ifelse(promoter.pqs$pqs.type=="pqs",ifelse(promoter.pqs$s2=="2","eG4","non-eG4"),"other")
#把tpm匹配到promoter.pqs，tpm17494,promoter17754,匹配不到的表达量设为0
promoter.pqs$kc.tpm <- tpm[match(promoter.pqs$id,tpm$gene_id),5]
#promoter.pqs[is.na(promoter.pqs$kc.tpm), 12]=0
promoter.pqs$s2.tpm <- tpm[match(promoter.pqs$id,tpm$gene_id),6]
#promoter.pqs[is.na(promoter.pqs$s2.tpm), 13]=0
#另一种
promoter.pqs <- subset(promoter.pqs,s2.tpm!='NA')
#画图
promoter.pqs$kc.type <- factor(promoter.pqs$kc.type,levels = c("non-eG4","eG4","other"))
promoter.pqs$s2.type <- factor(promoter.pqs$s2.type,levels = c("non-eG4","eG4","other"))
my_comparisons = list(c("non-eG4","eG4"),c("eG4","other"))
ggplot(data = promoter.pqs,aes(x=kc.type,y=log2(kc.tpm+1),fill=kc.type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(12.5,14),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="Promoter genes (Kc167 cells)") +
  coord_cartesian(ylim = c(0,15)) 
#ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"eG4_TPM_in_Kc.pdf"),
#       device = "pdf",width = 3.5,height = 3)
ggplot(data = promoter.pqs,aes(x=s2.type,y=log2(s2.tpm+1),fill=s2.type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(12.5,13.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#4d4d4d","#D7A49D","#92c5de")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="Promoter genes (S2 cells)") +
  coord_cartesian(ylim = c(0,15)) 
  
# ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"eG4_TPM_in_S2.pdf"),
#       device = "pdf",width = 3.5,height = 3)

#### 统计含有eG4基因表达的比例 ####
promoter.pqs$kc.tpm.type <- ifelse(promoter.pqs$kc.tpm>=1,"express","nonexpress")
table(promoter.pqs$kc.tpm.type,promoter.pqs$kc.type)
df <- table(promoter.pqs$kc.tpm.type,promoter.pqs$kc.type) %>% as.data.frame()
df <- spread(df,Var1,Freq)
df$sum <- rowSums(df[,2:3])
df$express_ratio <- df$express/df$sum

promoter.pqs$s2.tpm.type <- ifelse(promoter.pqs$s2.tpm>=1,"express","nonexpress")
table(promoter.pqs$s2.tpm.type,promoter.pqs$s2.type)
df <- table(promoter.pqs$s2.tpm.type,promoter.pqs$s2.type) %>% as.data.frame()
df <- spread(df,Var1,Freq)
df$sum <- rowSums(df[,2:3])
df$express_ratio <- df$express/df$sum

rm(list = ls());gc();rm(list = ls())#清空
Num = "007.4."
#### 分成两类not eG4,eG4, ####
promoter.kc.s2 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.2.promoter.kc.s2.txt") %>% as.data.frame()
promoter.kc.s2$kc.eG4 <- ifelse(promoter.kc.s2$kc==1,"eG4","no eG4")
promoter.kc.s2$s2.eG4 <- ifelse(promoter.kc.s2$s2==2,"eG4","no eG4")
tpm <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.Tpm.txt") %>% as.data.frame()
promoter.kc.s2$kc.tpm <- tpm[match(promoter.kc.s2$gene_id,tpm$gene_id),5]
promoter.kc.s2$s2.tpm <- tpm[match(promoter.kc.s2$gene_id,tpm$gene_id),6]
#promoter.kc.s2[is.na(promoter.kc.s2$kc.tpm),11]=0
#promoter.kc.s2[is.na(promoter.kc.s2$s2.tpm),12]=0
#另一种
promoter.kc.s2 <- subset(promoter.kc.s2,s2.tpm!='NA')
promoter.kc.s2$kc.eG4 <- factor(promoter.kc.s2$kc.eG4,levels = c("no eG4","eG4"))
promoter.kc.s2$s2.eG4 <- factor(promoter.kc.s2$s2.eG4,levels = c("no eG4","eG4"))
my_comparisons = list(c("no eG4","eG4"))
ggplot(data = promoter.kc.s2,aes(x=kc.eG4,y=log2(kc.tpm+1),fill=kc.eG4)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = 13.5,
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#4d4d4d","#FDDBC7")) +
  theme_bw()+
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="Promoter genes (Kc167 cells)") +
  coord_cartesian(ylim = c(0,15)) 
#ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"eG4_TPM_in_Kc.pdf"),
       device = "pdf",width = 3.5,height = 3)

ggplot(data = promoter.kc.s2,aes(x=s2.eG4,y=log2(s2.tpm+1),fill=s2.eG4)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = 13.5,
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+ ##默认使用Wilcoxon秩和检验
  scale_fill_manual(values = c("#4d4d4d","#D7A49D")) +
  theme_bw()+
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="Promoter genes (S2 cells)") +
  coord_cartesian(ylim = c(0,15)) 
#ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"eG4_TPM_in_Kc.pdf"),
device = "pdf",width = 3.5,height = 3)

#### 统计含有eG4基因表达的比例 ####
promoter.kc.s2$kc.tpm.type <- ifelse(promoter.kc.s2$kc.tpm>=1,"express","nonexpress")
table(promoter.kc.s2$kc.tpm.type,promoter.kc.s2$kc.eG4)
df <- table(promoter.kc.s2$kc.tpm.type,promoter.kc.s2$kc.eG4) %>% as.data.frame()
df <- spread(df,Var1,Freq)
df$sum <- rowSums(df[,2:3])
df$express_ratio <- df$express/df$sum

promoter.kc.s2$s2.tpm.type <- ifelse(promoter.kc.s2$s2.tpm>=1,"express","nonexpress")
table(promoter.kc.s2$s2.tpm.type,promoter.kc.s2$s2.eG4)
df <- table(promoter.kc.s2$s2.tpm.type,promoter.kc.s2$s2.eG4) %>% as.data.frame()
df <- spread(df,Var1,Freq)
df$sum <- rowSums(df[,2:3])
df$express_ratio <- df$express/df$sum

rm(list = ls());gc();rm(list = ls())#清空
Num = "007.4."
#*promoter区：基因上下游1000bp 基因含有eG4和不含有eG4的TPM---------------------------------------------
##1.promoter区：基因上下游1000bp
gene.sort <- fread('/home/yuss/flyG4/data/ref/gene.sort.bed') %>% as.data.frame()
gene.sort.p <- gene.sort[gene.sort$V5=='+',] ##提取正链的信息
gene.sort.m <- gene.sort[gene.sort$V5=='-',]
promoter.p <- data.frame(gene.sort.p$V1,gene.sort.p$V2-1000,gene.sort.p$V2+1000,gene.sort.p$V4,gene.sort.p$V5)
colnames(promoter.p) <- c("chr","start","end","gene_id","strand")
promoter.m <- data.frame(gene.sort.m$V1,gene.sort.m$V3-1000,gene.sort.m$V3+1000,gene.sort.m$V4,gene.sort.m$V5)
colnames(promoter.m) <- c("chr","start","end","gene_id","strand")
promoter <- rbind(promoter.p,promoter.m)
promoter$start <- ifelse(promoter$start <0, 0, promoter$start)
write.table(promoter,file = '/home/yuss/flyG4/data/ref/promoter1000.bed',
            sep = '\t',col.names = T,row.names = F,quote = F)
##2.promoter.bed与kc、s2细胞系中所有的G4取交集，若有交集表示启动子上的该基因片段有G4
# sed '1d' /home/yuss/flyG4/result/PQS/001.2.kc_all.bed > /home/yuss/flyG4/data/ref/kc_all.bed ##删除信息第一行
# sed '1d' /home/yuss/flyG4/result/PQS/001.2.s2_all.bed > /home/yuss/flyG4/data/ref/s2_all.bed
# sed '1d' /home/yuss/flyG4/data/ref/promoter1000.bed | bedtools intersect -a - -b /home/yuss/flyG4/data/ref/kc_all.bed -c > /home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter1000.kc.bed
# sed '1d' /home/yuss/flyG4/data/ref/promoter1000.bed | bedtools intersect -a - -b /home/yuss/flyG4/data/ref/s2_all.bed -c > /home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter1000.s2.bed

##3.读取promoter与s2,kc G4交集的文件
promoter.kc <- fread('/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter1000.kc.bed') %>% as.data.frame()
promoter.s2 <- fread('/home/yuss/flyG4/result//LucyCherbas.GR.2010.RNAseq/007.4.promoter1000.s2.bed') %>% as.data.frame()
#合并两个表
promoter.kc.s2 <- bind_cols(promoter.kc,promoter.s2$V6)
colnames(promoter.kc.s2) <- c("chr","start","end","gene_id","strand","kc","s2")
promoter.kc.s2$kc <- ifelse(promoter.kc.s2$kc==0,0,1)
promoter.kc.s2$s2 <- ifelse(promoter.kc.s2$s2==0,0,2)
#求sum
promoter.kc.s2$sum <- rowSums(promoter.kc.s2[,6:7])
write.table(promoter.kc.s2,file = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/",Num,"promoter1000.kc.s2.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)
# promoter.kcG4 <- promoter.kc.s2[promoter.kc.s2$sum==1,]
# promoter.s2G4 <- promoter.kc.s2[promoter.kc.s2$sum==2,]
# promoter.overlapG4 <- promoter.kc.s2[promoter.kc.s2$sum==3,]
# promoter.non_G4 <- promoter.kc.s2[promoter.kc.s2$sum==0,]
# promoter.mergeG4 <- promoter.kc.s2[promoter.kc.s2$sum!=0,]

#### 导入TPM ####
tpm <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.Tpm.txt") %>% as.data.frame()
promoter.kc.s2 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter1000.kc.s2.txt") %>% as.data.frame()
#把tpm匹配到promoter.kc.s2，tpm17494,promoter17754,匹配不到的表达量扔掉
promoter.kc.s2$kc.tpm <- tpm[match(promoter.kc.s2$gene_id,tpm$gene_id),5]
#promoter.pqs[is.na(promoter.pqs$kc.tpm), 12]=0
promoter.kc.s2$s2.tpm <- tpm[match(promoter.kc.s2$gene_id,tpm$gene_id),6]
#promoter.pqs[is.na(promoter.pqs$s2.tpm), 13]=0
#另一种
promoter.kc.s2 <- subset(promoter.kc.s2,s2.tpm!='NA')

# promoter.pqs$kc.tpm.type <- ifelse(promoter.pqs$kc.tpm>=1,"express","nonexpress")
# promoter.pqs$s2.tpm.type <- ifelse(promoter.pqs$s2.tpm>=1,"express","nonexpress")
# table(promoter.pqs$kc.tpm.type,promoter.pqs$kc.type)
#### 分成两类not eG4,eG4, ####
promoter.kc.s2$kc.eG4 <- ifelse(promoter.kc.s2$kc=="1","eG4","no eG4") #增加kctype列
promoter.kc.s2$s2.eG4 <- ifelse(promoter.kc.s2$s2=="2","eG4","no eG4") #增加kctype列
#画图
promoter.kc.s2$kc.eG4 <- factor(promoter.kc.s2$kc.eG4,levels = c("no eG4","eG4"))
promoter.kc.s2$s2.eG4 <- factor(promoter.kc.s2$s2.eG4,levels = c("no eG4","eG4"))
my_comparisons = list(c("no eG4","eG4"))
ggplot(data = promoter.kc.s2,aes(x=kc.eG4,y=log2(kc.tpm+1),fill=kc.eG4)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = 13.5,
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  theme_bw()+
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="Promoter genes (Kc167 cells)") +
  coord_cartesian(ylim = c(0,15)) 
#ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"eG4_TPM_in_Kc.pdf"),
device = "pdf",width = 3.5,height = 3)

ggplot(data = promoter.kc.s2,aes(x=s2.eG4,y=log2(s2.tpm+1),fill=s2.eG4)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = 13.5,
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+ ##默认使用Wilcoxon秩和检验
  scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  theme_bw()+
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="Promoter genes (S2 cells)") +
  coord_cartesian(ylim = c(0,15)) 
#ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"eG4_TPM_in_Kc.pdf"),
device = "pdf",width = 3.5,height = 3)


#00516C

#*promoter区：基因上下游2000bp 含有eG4和不含有eG4基因的TPM---------------------------------------------
rm(list = ls());gc();rm(list = ls())#清空
Num = "007.4."
##1.基因上下游2000bp
gene.sort <- fread('/home/yuss/flyG4/data/ref/gene.sort.bed') %>% as.data.frame()
gene.sort.p <- gene.sort[gene.sort$V5=='+',] ##提取正链的信息
gene.sort.m <- gene.sort[gene.sort$V5=='-',]
promoter.p <- data.frame(gene.sort.p$V1,gene.sort.p$V2-2000,gene.sort.p$V2+2000,gene.sort.p$V4,gene.sort.p$V5)
colnames(promoter.p) <- c("chr","start","end","gene_id","strand")
promoter.m <- data.frame(gene.sort.m$V1,gene.sort.m$V3-2000,gene.sort.m$V3+2000,gene.sort.m$V4,gene.sort.m$V5)
colnames(promoter.m) <- c("chr","start","end","gene_id","strand")
promoter <- rbind(promoter.p,promoter.m)
promoter$start <- ifelse(promoter$start <0, 0, promoter$start)
write.table(promoter,file = '/home/yuss/flyG4/data/ref/promoter2000.bed',
            sep = '\t',col.names = T,row.names = F,quote = F)
##2.promoter.bed与kc、s2细胞系中所有的G4取交集，若有交集表示启动子上的该基因片段有G4
# sed '1d' /home/yuss/flyG4/result/PQS/001.2.kc_all.bed > /home/yuss/flyG4/data/ref/kc_all.bed ##删除信息第一行
# sed '1d' /home/yuss/flyG4/result/PQS/001.2.s2_all.bed > /home/yuss/flyG4/data/ref/s2_all.bed
# sed '1d' /home/yuss/flyG4/data/ref/promoter2000.bed | bedtools intersect -a - -b /home/yuss/flyG4/data/ref/kc_all.bed -c > /home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter2000.kc.bed
# sed '1d' /home/yuss/flyG4/data/ref/promoter2000.bed | bedtools intersect -a - -b /home/yuss/flyG4/data/ref/s2_all.bed -c > /home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter2000.s2.bed

##3.读取promoter与s2,kc G4交集的文件
promoter.kc <- fread('/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter2000.kc.bed') %>% as.data.frame()
promoter.s2 <- fread('/home/yuss/flyG4/result//LucyCherbas.GR.2010.RNAseq/007.4.promoter2000.s2.bed') %>% as.data.frame()
#合并两个表
promoter.kc.s2 <- bind_cols(promoter.kc,promoter.s2$V6)
colnames(promoter.kc.s2) <- c("chr","start","end","gene_id","strand","kc","s2")
promoter.kc.s2$kc <- ifelse(promoter.kc.s2$kc==0,0,1)
promoter.kc.s2$s2 <- ifelse(promoter.kc.s2$s2==0,0,2)
#求sum
promoter.kc.s2$sum <- rowSums(promoter.kc.s2[,6:7])
write.table(promoter.kc.s2,file = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/",Num,"promoter2000.kc.s2.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)
# promoter.kcG4 <- promoter.kc.s2[promoter.kc.s2$sum==1,]
# promoter.s2G4 <- promoter.kc.s2[promoter.kc.s2$sum==2,]
# promoter.overlapG4 <- promoter.kc.s2[promoter.kc.s2$sum==3,]
# promoter.non_G4 <- promoter.kc.s2[promoter.kc.s2$sum==0,]
# promoter.mergeG4 <- promoter.kc.s2[promoter.kc.s2$sum!=0,]
# #按行合并

#### 导入TPM ####
tpm <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.Tpm.txt") %>% as.data.frame()
promoter.kc.s2 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter2000.kc.s2.txt") %>% as.data.frame()
#把tpm匹配到promoter.kc.s2，tpm17494,promoter17754,匹配不到的表达量扔掉
promoter.kc.s2$kc.tpm <- tpm[match(promoter.kc.s2$gene_id,tpm$gene_id),5]
#promoter.kc.s2[is.na(promoter.kc.s2$kc.tpm), 9]=0
promoter.kc.s2$s2.tpm <- tpm[match(promoter.kc.s2$gene_id,tpm$gene_id),6]
#promoter.kc.s2[is.na(promoter.kc.s2$s2.tpm), 10]=0
#另一种
promoter.kc.s2 <- subset(promoter.kc.s2,s2.tpm!='NA')

#### 分成两类no eG4,eG4 ####
promoter.kc.s2$kc.eG4 <- ifelse(promoter.kc.s2$kc=="1","eG4","no eG4") #增加kctype列
promoter.kc.s2$s2.eG4 <- ifelse(promoter.kc.s2$s2=="2","eG4","no eG4") #增加kctype列
#画图
promoter.kc.s2$kc.eG4 <- factor(promoter.kc.s2$kc.eG4,levels = c("no eG4","eG4"))
promoter.kc.s2$s2.eG4 <- factor(promoter.kc.s2$s2.eG4,levels = c("no eG4","eG4"))
my_comparisons = list(c("no eG4","eG4"))
ggplot(data = promoter.kc.s2,aes(x=kc.eG4,y=log2(kc.tpm+1),fill=kc.eG4)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     method="t.test",
                     paired = FALSE,  # 设置为 FALSE，表示独立样本 t 检验
                     label.y = 13.5,
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#8EA1CC","#FC8E62")) +
  theme_bw()+
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2]~TPM))+
  labs(title="Promoter genes (Kc167 cells)") +
  coord_cartesian(ylim = c(0,15)) 
ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"Promoter.noeG4.eG4_TPM_in_Kc.pdf"),
device = "pdf",width = 3.3,height = 4.2)

ggplot(data = promoter.kc.s2,aes(x=s2.eG4,y=log2(s2.tpm+1),fill=s2.eG4)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     method="t.test",
                     paired = FALSE,  # 设置为 FALSE，表示独立样本 t 检验
                     label.y = 13.5,
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+ ##默认使用Wilcoxon秩和检验
  scale_fill_manual(values = c("#8EA1CC","#66C3AA")) +
  theme_bw()+
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2]~TPM))+
  labs(title="Promoter genes (S2 cells)") +
  coord_cartesian(ylim = c(0,15)) 
ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"Promoter.noeG4.eG4_TPM_in_S2.pdf"),
       device = "pdf",width = 3.3,height = 4.2)
aggregate(s2.tpm~ s2.eG4, data=promoter.kc.s2, FUN=median)
pro.s2.eG4 <- promoter.kc.s2 %>% filter(s2.eG4=="eG4")
pro.s2.noeG4 <- promoter.kc.s2 %>% filter(s2.eG4=="no eG4")
t.test(pro.s2.noeG4$s2.tpm,pro.s2.eG4$s2.tpm)
# Welch Two Sample t-test
# 
# data:  pro.s2.noeG4$s2.tpm and pro.s2.eG4$s2.tpm
# t = 0.60669, df = 6073.4, p-value = 0.5441
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -13.84876  26.26244
# sample estimates:
#   mean of x mean of y 
# 76.05825  69.85141 

#### promoter区所有染色体两类G4的表达####
my_comparisons <- list(c("no eG4", "eG4"))
subset_promoter.kc.s2 <- subset(promoter.kc.s2, chr %in% c("2L", "2R", "3L","3R","X"))
ggplot(data = subset_promoter.kc.s2,aes(x=kc.eG4,y=log2(kc.tpm+1),fill=kc.eG4)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  facet_grid(~ chr) +
  stat_compare_means(comparisons = my_comparisons,
                     method="t.test",
                     paired = FALSE,  # 设置为 FALSE，表示独立样本 t 检验
                     label.y = c(rep(13,5)),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#8EA1CC","#FC8E62")) +
  theme_bw()+
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2]~TPM))+
  labs(title="Promoter genes (Kc167 cells)") +
  coord_cartesian(ylim = c(0,15)) #

ggplot(data = subset_promoter.kc.s2,aes(x=s2.eG4,y=log2(s2.tpm+1),fill=s2.eG4)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  facet_grid(~ chr) +
  stat_compare_means(comparisons = my_comparisons,
                     method="t.test",
                     paired = FALSE,  # 设置为 FALSE，表示独立样本 t 检验
                     label.y = c(rep(13,5)),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#8EA1CC","#66C3AA")) +
  theme_bw()+
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2]~TPM))+
  labs(title="Promoter genes (S2 cells)") +
  coord_cartesian(ylim = c(0,15))   
aggregate(s2.tpm ~ s2.eG4, data = subset_promoter.kc.s2, FUN = mean)

result <- subset_promoter.kc.s2 %>%
  filter(chr == "X") %>%
  group_by(s2.eG4) %>%
  summarize(mean_s2.tpm = mean(s2.tpm))


#### promoter区X染色体上两类的表达水平 ####
promoter.kc.s2.X <- promoter.kc.s2[promoter.kc.s2$chr=="X",]
ggplot(data = promoter.kc.s2.X,aes(x=kc.eG4,y=log2(kc.tpm+1),fill=kc.eG4)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = 13.5,
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#8EA1CC","#FC8E62")) +
  theme_bw()+
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2]~TPM))+
  labs(title="Promoter genes (Kc167 cells)") +
  coord_cartesian(ylim = c(0,15)) 
#ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"Promoter.noeG4.eG4_TPM_in_Kc.pdf"),
       device = "pdf",width = 3.3,height = 4.2)

ggplot(data = promoter.kc.s2.X,aes(x=s2.eG4,y=log2(s2.tpm+1),fill=s2.eG4)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = 13.5,
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+ ##默认使用Wilcoxon秩和检验
  scale_fill_manual(values = c("#8EA1CC","#66C3AA")) +
  theme_bw()+
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2]~TPM))+
  labs(title="Promoter genes (S2 cells)") +
  coord_cartesian(ylim = c(0,15)) 
#ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"Promoter.noeG4.eG4_TPM_in_S2.pdf"),
       device = "pdf",width = 3.3,height = 4.2)

#### X染色体上含有G4的基因的表达比例 ####
promoter.kc.s2.X$kc.tpm.type <- ifelse(promoter.kc.s2.X$kc.tpm>=1,"express","nonexpress")

# X.gene.kc.s2$a <- promoter.kc.s2.X[match(X.gene.kc.s2$id,promoter.kc.s2.X$gene_id),1]
# # 检查某一列是否有NA值
# if (any(is.na(X.gene.kc.s2$a))) {
#   cat("Column 'a' contains NA values.\n")
# } else {
#   cat("Column 'a' does not contain any NA values.\n")
# }
# 所有X染色体上的基因都位于promoter区

table(promoter.kc.s2.X$kc.tpm.type,promoter.kc.s2.X$kc.eG4)
kc.df <- table(promoter.kc.s2.X$kc.tpm.type,promoter.kc.s2.X$kc.eG4) %>% as.data.frame()
#长数据转为宽数据
kc.df <- spread(kc.df,Var1,Freq)
kc.df$sum <- rowSums(kc.df[,2:3])
kc.df$express_ratio <- kc.df$express/kc.df$sum
kc.s2.express <- kc.df[,c(1,5)]
promoter.kc.s2.X$s2.tpm.type <- ifelse(promoter.kc.s2.X$s2.tpm>=1,"express","nonexpress")
table(promoter.kc.s2.X$s2.tpm.type,promoter.kc.s2.X$s2.eG4)
s2.df <- table(promoter.kc.s2.X$s2.tpm.type,promoter.kc.s2.X$s2.eG4) %>% as.data.frame()
#长数据转为宽数据
s2.df <- spread(s2.df,Var1,Freq)
s2.df$sum <- rowSums(s2.df[,2:3])
s2.df$express_ratio <- s2.df$express/s2.df$sum
kc.s2.express$s2 <- s2.df[,5]
colnames(kc.s2.express)[2] <- "kc"
kc.s2.express.long <- pivot_longer(kc.s2.express, cols = c(kc, s2), names_to = "Measure", values_to = "Value")
kc.s2.express.long <- kc.s2.express.long %>% mutate(type = paste(Var2, Measure,sep = "_"))
kc.s2.express.long$type <- factor(kc.s2.express.long$type,levels = c("no eG4_kc","eG4_kc","no eG4_s2","eG4_s2"))
kc.s2.express.long$Value <- round(kc.s2.express.long$Value,2)*100
kc.s2.express.long$percent <- paste(kc.s2.express.long$Value, "%",sep = "")


ggplot(kc.s2.express.long,aes(x=type,y=Value,fill=type))+
  geom_bar(stat = "identity",position = "dodge",width = 0.8)+
  theme_bw()+ 
  ylab("% Expressed genes")+
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=10,colour = "black"),axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,colour = "black"),axis.title.y = element_text(size=12),
        legend.position = "none")+
  scale_y_continuous(expand = c(0,0),limits = c(0,70)) +
  geom_vline(xintercept = 2.5,color="#595959",size=0.4) +
  geom_text(aes(label=percent,vjust=-0.5),color="black",size=4) + 
  labs(title="Promoter genes on the X chromosome (Kc167 and S2 cells)") +
  scale_fill_manual(values = c("no eG4_kc"="#8EA1CC","eG4_kc"="#FC8E62","no eG4_s2"="#8EA1CC","eG4_s2"="#66C3AA")) +
  scale_x_discrete(labels = c("no eG4","eG4","no eG4","eG4"))
table(promoter.kc.s2.X$kc.tpm.type)
table(X.gene.kc.s2$kc.tpm.type)
table(promoter.kc.s2.X$s2.tpm.type)
table(X.gene.kc.s2$s2.tpm.type)
#### 统计含有eG4基因表达的比例 ####
promoter.kc.s2$kc.tpm.type <- ifelse(promoter.kc.s2$kc.tpm>=1,"express","nonexpress")
table(promoter.kc.s2$kc.tpm.type,promoter.kc.s2$kc.eG4)
df <- table(promoter.kc.s2$kc.tpm.type,promoter.kc.s2$kc.eG4) %>% as.data.frame()
df <- spread(df,Var1,Freq)
df$sum <- rowSums(df[,2:3])
df$express_ratio <- df$express/df$sum
kc.s2.express <- df[,c(1,5)]
promoter.kc.s2$s2.tpm.type <- ifelse(promoter.kc.s2$s2.tpm>=1,"express","nonexpress")
table(promoter.kc.s2$s2.tpm.type,promoter.kc.s2$s2.eG4)
df <- table(promoter.kc.s2$s2.tpm.type,promoter.kc.s2$s2.eG4) %>% as.data.frame()
df <- spread(df,Var1,Freq)
df$sum <- rowSums(df[,2:3])
df$express_ratio <- df$express/df$sum
kc.s2.express$s2 <- df[,5]
colnames(kc.s2.express)[2] <- "kc"
kc.s2.express.long <- pivot_longer(kc.s2.express, cols = c(kc, s2), names_to = "Measure", values_to = "Value")
kc.s2.express.long <- kc.s2.express.long %>% mutate(type = paste(Var2, Measure,sep = "_"))
kc.s2.express.long$type <- factor(kc.s2.express.long$type,levels = c("no eG4_kc","eG4_kc","no eG4_s2","eG4_s2"))
kc.s2.express.long$Value <- round(kc.s2.express.long$Value,2)*100
kc.s2.express.long$percent <- paste(kc.s2.express.long$Value, "%",sep = "")
ggplot(kc.s2.express.long,aes(x=type,y=Value,fill=type))+
  geom_bar(stat = "identity",position = "dodge",width = 0.8)+
  theme_bw()+ 
  ylab("% Expressed genes")+
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=10,colour = "black"),axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,colour = "black"),axis.title.y = element_text(size=12),
        legend.position = "none")+
  scale_y_continuous(expand = c(0,0),limits = c(0,70)) +
  geom_vline(xintercept = 2.5,color="#595959",size=0.4) +
  geom_text(aes(label=percent,vjust=-0.5),color="black",size=4) + 
  labs(title="Promoter genes (Kc167 and S2 cells)") +
  scale_fill_manual(values = c("no eG4_kc"="#8EA1CC","eG4_kc"="#FC8E62","no eG4_s2"="#8EA1CC","eG4_s2"="#66C3AA")) +
  scale_x_discrete(labels = c("no eG4","eG4","no eG4","eG4"))
ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"Promoter.noeG4.eG4.ExpressRatio_in_KcS2.pdf"),
       device = "pdf",width = 4,height = 4.8)

rm(list = ls());gc();rm(list = ls())#清空
Num = "007.4."

#### 分成三类other,eG4,non-eG4 ####
promoter.pqs <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter2000.pqs.bed") %>% as.data.frame()
colnames(promoter.pqs) <- c("chr","start","end","id","strand","promoter.pqs")
promoter.pqs$pqs.type <- ifelse(promoter.pqs$promoter.pqs==0,"other","pqs")
promoter.kc.s2 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.2.promoter.kc.s2.txt") %>% as.data.frame()
promoter.pqs$kc <- promoter.kc.s2[match(promoter.pqs$id,promoter.kc.s2$gene_id),6] #匹配kc列
promoter.pqs$kc.type <- ifelse(promoter.pqs$pqs.type=="pqs",ifelse(promoter.pqs$kc=="1","eG4","non-eG4"),"other") #增加kctype列
promoter.pqs$s2 <- promoter.kc.s2[match(promoter.pqs$id,promoter.kc.s2$gene_id),7] #匹配s2列
promoter.pqs$s2.type <- ifelse(promoter.pqs$pqs.type=="pqs",ifelse(promoter.pqs$s2=="2","eG4","non-eG4"),"other")
#把tpm匹配到promoter.pqs，tpm17494,promoter17754,匹配不到的表达量设为0
tpm <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.Tpm.txt") %>% as.data.frame()
promoter.pqs$kc.tpm <- tpm[match(promoter.pqs$id,tpm$gene_id),5]
#promoter.pqs[is.na(promoter.pqs$kc.tpm), 12]=0
promoter.pqs$s2.tpm <- tpm[match(promoter.pqs$id,tpm$gene_id),6]
#promoter.pqs[is.na(promoter.pqs$s2.tpm), 13]=0
#另一种
promoter.pqs <- subset(promoter.pqs,s2.tpm!='NA')
#画图
promoter.pqs$kc.type <- factor(promoter.pqs$kc.type,levels = c("non-eG4","eG4","other"))
promoter.pqs$s2.type <- factor(promoter.pqs$s2.type,levels = c("non-eG4","eG4","other"))
my_comparisons = list(c("non-eG4","eG4"),c("eG4","other"))
ggplot(data = promoter.pqs,aes(x=kc.type,y=log2(kc.tpm+1),fill=kc.type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(12.5,14),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="Promoter genes (Kc167 cells)") +
  coord_cartesian(ylim = c(0,15)) 
#ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"eG4_TPM_in_Kc.pdf"),
#       device = "pdf",width = 3.5,height = 3)
ggplot(data = promoter.pqs,aes(x=s2.type,y=log2(s2.tpm+1),fill=s2.type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(12.5,13.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#4d4d4d","#D7A49D","#92c5de")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="Promoter genes (S2 cells)") +
  coord_cartesian(ylim = c(0,15)) 

# ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"eG4_TPM_in_S2.pdf"),
#       device = "pdf",width = 3.5,height = 3)

#### 统计含有eG4基因表达的比例 ####
promoter.pqs$kc.tpm.type <- ifelse(promoter.pqs$kc.tpm>=1,"express","nonexpress")
table(promoter.pqs$kc.tpm.type,promoter.pqs$kc.type)
df <- table(promoter.pqs$kc.tpm.type,promoter.pqs$kc.type) %>% as.data.frame()
df <- spread(df,Var1,Freq)
df$sum <- rowSums(df[,2:3])
df$express_ratio <- df$express/df$sum

promoter.pqs$s2.tpm.type <- ifelse(promoter.pqs$s2.tpm>=1,"express","nonexpress")
table(promoter.pqs$s2.tpm.type,promoter.pqs$s2.type)
df <- table(promoter.pqs$s2.tpm.type,promoter.pqs$s2.type) %>% as.data.frame()
df <- spread(df,Var1,Freq)
df$sum <- rowSums(df[,2:3])
df$express_ratio <- df$express/df$sum


my_comparisons <- list(c("non-eG4", "eG4"),c("eG4","other"))
subset_promoter.pqs <- subset(promoter.pqs, chr %in% c("2L", "2R", "3L","3R","X"))
ggplot(data = subset_promoter.pqs,aes(x=kc.type,y=log2(kc.tpm+1),fill=kc.type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  facet_grid(~ chr) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = rep(c(12.5,13.5),5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  theme_bw()+
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2]~TPM))+
  labs(title="Promoter genes (Kc167 cells)") +
  coord_cartesian(ylim = c(0,15)) #


ggplot(data = subset_promoter.pqs,aes(x=s2.type,y=log2(s2.tpm+1),fill=s2.type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  facet_grid(~ chr) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = rep(c(12.5,13.5),5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#4d4d4d","#D7A49D","#92c5de")) +
  theme_bw()+
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2]~TPM))+
  labs(title="Promoter genes (Kc167 cells)") +
  coord_cartesian(ylim = c(0,15)) #


aggregate(s2.tpm ~ s2.eG4, data = subset_promoter.kc.s2, FUN = median)