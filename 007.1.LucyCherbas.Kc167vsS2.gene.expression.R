rm(list = ls());gc();rm(list = ls())#清空
Num = "007.1."
getwd()
setwd("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/")
#### counts ####
summary_readcounts <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/merged_counts.txt") %>% as.data.frame()
readcounts <- summary_readcounts[7:ncol(summary_readcounts)]
colnames(readcounts) <- unlist(lapply(strsplit(colnames(readcounts), "[/]"),"[[",8))
# strsplit(colnames(readcounts), "[/]") 使用 strsplit 函数，根据斜杠字符 ("/") 将每个列名分割成一个字符串向量。这将生成一个列表，列表中的每个元素是一个字符串向量，包含分割后的各部分。
# lapply(..., "[[", 8) 使用 lapply 函数来遍历上一步生成的列表，并提取每个字符串向量的第 8 个元素，这是分割后的文件名部分。
# unlist(...) 将前一步得到的提取文件名的结果合并成一个字符向量。
colnames(readcounts) <- unlist(lapply(strsplit(colnames(readcounts),"[.]"),"[[",1))
row.names(readcounts) <- summary_readcounts$Geneid
readcounts<-readcounts[apply(readcounts, 1, sum)>0,]
counts <- readcounts
counts$geneid <- rownames(counts)
counts <- counts[,c(5,1,2,3,4)]
write.table(counts,file = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/",Num,"counts.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)
sample_name <- factor(colnames(readcounts))
metadata <- data.frame(sample_name)
metadata$group <- as.factor(c("kc","kc","s2","s2"))

#### Preparation for analysis ####
library(DESeq2)
##1.构建矩阵
dds <- DESeqDataSetFromMatrix(countData = readcounts,
                              colData = metadata,
                              design = ~ group)
##2.dds标准化
dds <- DESeq(dds)
##3.获取标准化后的数据
normalized_counts <- counts(dds,normalized=T)
##4.有两个重复取平均
normalized_mean_counts = t(apply(normalized_counts, 1, function(a){tapply(a, metadata$group, mean)}))
write.table(normalized_counts,file = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/",Num,"normalized_counts.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)
write.table(normalized_mean_counts,file = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/",Num,"normalized_mean_counts.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)
#### 差异表达分析 ####
res <- results(dds,contrast = c("group","kc","s2")) ##谁比谁
summary(res)
res1 <- data.frame(res,stringsAsFactors = F,check.names = F) %>% na.omit()
res1 <- cbind(res1, normalized_mean_counts[rownames(res1), c("kc", "s2")])
res1$signifi <- ifelse(res1$pvalue<0.05&abs(res1$log2FoldChange)>=0,ifelse(res1$log2FoldChange>0,"Up","Down"),"No-sig")
table(res1$signifi)
res1 <- rownames_to_column(res1, 'gene_id') 
write.table(res1,file = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/",Num,"DE.kcvss2.txt"),
sep = '\t',col.names = T,row.names = F,quote = F)
#5.火山图可视化
library(ggplot2)
ggplot(res1,aes(x=log2FoldChange,y=-log10(pvalue))) +#黑白主题，base_size可以设置某种主题内基本字体的大小
  geom_point(aes(color=signifi),alpha=0.5,size=2) +
  theme_bw(base_size = 16)+
  theme(#aspect.ratio = 1,
    #plot.title = element_text(hjust = 0.5), # #设置标题居中
    axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12),
    panel.grid.major = element_blank(),  # 去除主要网格
    panel.grid.minor = element_blank(),  # 去除次要网格
    panel.background = element_rect(fill = "white"))+
  scale_color_manual(name = '',
                     values = c('Up'='#D6604D','Nosig'='grey','Down'='#74ADD1'), #手动设置颜色时调整颜色的因子顺序
                     label = c('Up'='KC','NoSig'='Nosig','Down'='S2'))+
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 0.8) +
  geom_vline(xintercept = c(-1,1),lty = 'dashed',size = 0.8) +
  annotate(geom="text", x=-6.5, y=270, size=4,label=paste0("S2\n(N=",sum(res1$signifi=="Down"),")"))+
  annotate(geom="text", x=6.5, y=270, size=4,label=paste0("Kc167\n(N=",sum(res1$signifi=="Up"),")"))+
  guides(color="none") + #图例颜色
  #scale_x_continuous(breaks = c(-6,-4,-2,-1,0,1,2,4,6)) +
  ## 修改坐标轴
  xlab(bquote(Log[2]~Fold~Change))+
  ylab(bquote(Log[10]~P~Value))
ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"DEKc167vsS2.VolcanoPlot.pdf"),
       device = "pdf",width = 4,height = 3) 

#### 差异基因和非差异基因mergeG4的富集####
tpm <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.Tpm.txt") %>% as.data.frame()
res1$kc.tpm <- tpm[match(res1$gene_id,tpm$gene_id),5]
res1$s2.tpm <- tpm[match(res1$gene_id,tpm$gene_id),6]

gene <- fread("/home/yuss/flyG4/data/ref/gene.bed") %>% as.data.frame()
res1$chr <- gene[match(res1$gene_id,gene$V4),1]
res1$start <- gene[match(res1$gene_id,gene$V4),2]
res1$end <- gene[match(res1$gene_id,gene$V4),3]
up <- res1[res1$signifi=="Up",]
down <- res1[res1$signifi=="Down",]
DEG <- res1[res1$signifi!="No-sig",]
nonDEG <- res1[res1$signifi=="No-sig",]
DEG <- DEG[,c(11,14:15,1)]
nonDEG <- nonDEG[,c(11,14:15,1)]
write.table(DEG,file = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/",Num,"DEG.bed"),
            sep = '\t',col.names = F,row.names = F,quote = F)
write.table(nonDEG,file = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/",Num,"nonDEG.bed"),
            sep = '\t',col.names = F,row.names = F,quote = F)

DEG <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.1.DEG.mergeG4.bed") %>% as.data.frame()
sum(DEG$V5=="0")##DEG没有G4的数量3401
sum(DEG$V5!="0")##DEG有G4的数量1455
nonDEG <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.1.nonDEG.mergeG4.bed") %>% as.data.frame()
sum(nonDEG$V5=="0")##非DEG没有G4的数量4094
sum(nonDEG$V5!="0")##非DEG有G4的数量1334
#### fishertest ####
data <- matrix(c(1445,1334,3401,4094),nrow = 2)
colnames(data) <- c("eG4","no eG4")
rownames(data) <- c("DEG","nonDEG")
fisher.test(data) # p-value = 2.459e-09
#看到了在差异基因上G4的富集

#### 上下调差异表达基因对应细胞系所有eG4的富集 ####
up <- up[,c(11,14:15,1)]
down <- down[,c(11,14:15,1)]
write.table(up,file = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/",Num,"upDEG.bed"),
            sep = '\t',col.names = F,row.names = F,quote = F)
write.table(down,file = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/",Num,"downDEG.bed"),
            sep = '\t',col.names = F,row.names = F,quote = F)

##交集
# (base) yuss@ubuntu:~/flyG4/result/LucyCherbas.GR.2010.RNAseq$ bedtools intersect -a 007.1.upDEG.bed -b "/home/yuss/flyG4/result/PQS/5.type.bed/001.2.kc_all.bed" -c > 007.1.upDEG.kcallG4.bed
# (base) yuss@ubuntu:~/flyG4/result/LucyCherbas.GR.2010.RNAseq$ bedtools intersect -a 007.1.upDEG.bed -b "/home/yuss/flyG4/result/PQS/5.type.bed/001.2.kc_specific.bed" -c > 007.1.upDEG.kcspecificG4.bed
# (base) yuss@ubuntu:~/flyG4/result/LucyCherbas.GR.2010.RNAseq$ bedtools intersect -a 007.1.downDEG.bed -b "/home/yuss/flyG4/result/PQS/5.type.bed/001.2.s2_all.bed" -c > 007.1.downDEG.s2allG4.bed
# (base) yuss@ubuntu:~/flyG4/result/LucyCherbas.GR.2010.RNAseq$ bedtools intersect -a 007.1.downDEG.bed -b "/home/yuss/flyG4/result/PQS/5.type.bed/001.2.s2_specific.bed" -c > 007.1.downDEG.s2specificG4.bed
# (base) yuss@ubuntu:~/flyG4/result/LucyCherbas.GR.2010.RNAseq$ bedtools intersect -a 007.1.nonDEG.bed -b "/home/yuss/flyG4/result/PQS/5.type.bed/001.2.kc_all.bed" -c > 007.1.nonDEG.kcallG4.bed
# (base) yuss@ubuntu:~/flyG4/result/LucyCherbas.GR.2010.RNAseq$ bedtools intersect -a 007.1.nonDEG.bed -b "/home/yuss/flyG4/result/PQS/5.type.bed/001.2.kc_specific.bed" -c > 007.1.nonDEG.kcspecificG4.bed
# (base) yuss@ubuntu:~/flyG4/result/LucyCherbas.GR.2010.RNAseq$ bedtools intersect -a 007.1.nonDEG.bed -b "/home/yuss/flyG4/result/PQS/5.type.bed/001.2.s2_all.bed" -c >  007.1.nonDEG.s2allG4.bed
# (base) yuss@ubuntu:~/flyG4/result/LucyCherbas.GR.2010.RNAseq$ bedtools intersect -a 007.1.nonDEG.bed -b "/home/yuss/flyG4/result/PQS/5.type.bed/001.2.s2_specific.bed" -c > 007.1.nonDEG.s2specificG4.bed

upG4 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.1.upDEG.kcallG4.bed") %>% as.data.frame()
nonDEG.kcallG4 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.1.nonDEG.kcallG4.bed") %>% as.data.frame()
sum(upG4$V5!=0) #上调表达的差异基因上有kcG4的数量746
sum(upG4$V5==0) #上调表达的差异基因上没有kcG4的数量1670
sum(nonDEG.kcallG4$V5!=0) #1219
sum(nonDEG.kcallG4$V5==0) #4209

data <- matrix(c(746,1219,1670,4209),nrow = 2)
colnames(data) <- c("eG4","no eG4")
rownames(data) <- c("DEG","nonDEG")
fisher.test(data) # p-value = 4.755e-15
#看到了在上调表达的差异基因上kcG4的富集
# > 746/(746+1670)
# [1] 0.3087748
# > 1219/(1219+4209)
# [1] 0.2245763

#### 上下调差异表达基因对应细胞系特有eG4的富集 ####
downG4 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.1.downDEG.s2allG4.bed") %>% as.data.frame()
nonDEG.s2allG4 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.1.nonDEG.s2allG4.bed") %>% as.data.frame()
sum(downG4$V5!=0) #下调表达的差异基因上有s2G4的数量558
sum(downG4$V5==0) #下调表达的差异基因上没有s2G4的数量1882
sum(nonDEG.s2allG4$V5!=0) #1133
sum(nonDEG.s2allG4$V5==0) #4295

data <- matrix(c(558,1133,1882,4295),nrow = 2)
colnames(data) <- c("eG4","no eG4")
rownames(data) <- c("DEG","nonDEG")
fisher.test(data) # p-value = 0.04688
#看到了在下调表达的差异基因上s2G4的富集
# > 558/(558+1133)
# [1] 0.3299823
# > 1882/(1882+4295)
# [1] 0.3046786

upG4 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.1.upDEG.kcspecificG4.bed") %>% as.data.frame()
nonDEG.kcspecificG4 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.1.nonDEG.kcspecificG4.bed") %>% as.data.frame()
sum(upG4$V5!=0) #上调表达的差异基因上有kcG4的数量323
sum(upG4$V5==0) #上调表达的差异基因上没有kcG4的数量2093
sum(nonDEG.kcspecificG4$V5!=0) #483
sum(nonDEG.kcspecificG4$V5==0) #4945

data <- matrix(c(323,483,2093,4945),nrow = 2)
colnames(data) <- c("eG4","no eG4")
rownames(data) <- c("DEG","nonDEG")
fisher.test(data) # p-value = 3.694e-09
#看到了在上调表达的差异基因上kcspecificG4的富集

downG4 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.1.downDEG.s2specificG4.bed") %>% as.data.frame()
nonDEG.s2specificG4 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.1.nonDEG.s2specificG4.bed") %>% as.data.frame()
sum(downG4$V5!=0) #下调表达的差异基因上有s2G4的数量174
sum(downG4$V5==0) #下调表达的差异基因上没有s2G4的数量2266
sum(nonDEG.s2specificG4$V5!=0) #355
sum(nonDEG.s2specificG4$V5==0) #5073

data <- matrix(c(174,355,2266,5073),nrow = 2)
colnames(data) <- c("eG4","no eG4")
rownames(data) <- c("DEG","nonDEG")
fisher.test(data) # p-value = 0.3308
#看到了在下调表达的差异基因上s2G4的富集

##用specific eG4去做，因为相同的eG4有差异表达基因，岂不是说明了差异表达不是eG4影响的
##应该kc\s2所有的eG4做，因为如果只用specific eG4做的话，那就认定了在不同细胞系的相同位置形成的eG4对基因表达调控是相同的
##但是不同细胞系的eG4对基因表达调控的作用是否一样，是不知道的