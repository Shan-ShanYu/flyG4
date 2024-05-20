rm(list = ls());gc();rm(list = ls())#清空
Num = "010.1."
#### counts ####
summary_readcounts <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/merged_counts.txt") %>% as.data.frame()
readcounts <- summary_readcounts[7:ncol(summary_readcounts)]
colnames(readcounts) <- unlist(lapply(strsplit(colnames(readcounts), "[/]"),"[[",8))
# strsplit(colnames(readcounts), "[/]") 使用 strsplit 函数，根据斜杠字符 ("/") 将每个列名分割成一个字符串向量。这将生成一个列表，列表中的每个元素是一个字符串向量，包含分割后的各部分。
# lapply(..., "[[", 8) 使用 lapply 函数来遍历上一步生成的列表，并提取每个字符串向量的第 8 个元素，这是分割后的文件名部分。
# unlist(...) 将前一步得到的提取文件名的结果合并成一个字符向量。
colnames(readcounts) <- unlist(lapply(strsplit(colnames(readcounts),"[.]"),"[[",1))
colnames(readcounts) <- sub("_1_val_1", "", colnames(readcounts))
row.names(readcounts) <- summary_readcounts$Geneid
readcounts<-readcounts[apply(readcounts, 1, sum)>0,]
readcounts <- readcounts[,c(1:8,10:17)]
sample_name <- factor(colnames(readcounts))
metadata <- data.frame(sample_name)
metadata$treat <- as.factor(rep(rep(c("con", "PDS", "Phen"), c(3, 3, 2)),times=2))
metadata$cell <- as.factor(rep(c("kc","s2"),each=8))
readcounts$geneid <- rownames(readcounts)

##分细胞系
readcounts_kc <- readcounts[,c(17,1:8)]
readcounts_s2 <- readcounts[,c(17,9:16)]
metadata_kc <- metadata[1:8,]
metadata_s2 <- metadata[9:16,]

write.table(readcounts,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"counts.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)
write.table(metadata,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"metadata.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)
write.table(readcounts_kc,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"kc.counts.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)
write.table(readcounts_s2,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"s2.counts.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)
write.table(metadata_kc,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"metadata_kc.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)
write.table(metadata_s2,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"metadata_s2.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)

#### Kc Phen\PDS vs con DESeq2 #### 
rm(list = ls());gc();rm(list = ls())#清空
Num = "010.1."
readcounts_kc <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.kc.counts.txt") %>% as.data.frame()
rownames(readcounts_kc) <- readcounts_kc$geneid
readcounts_kc <- readcounts_kc[,-1]
metadata_kc <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.metadata_kc.txt") %>% as.data.frame()
metadata_kc$treat <- as.factor(metadata_kc$treat)
library(DESeq2)
##1.构建矩阵
dds <- DESeqDataSetFromMatrix(countData = readcounts_kc,
                              colData = metadata_kc,
                              design = ~ treat)
##2.dds标准化
dds <- DESeq(dds)
##3.获取标准化后的数据
normalized_counts <- counts(dds,normalized=T)
##4.有三个重复取平均
normalized_mean_counts = t(apply(normalized_counts, 1, function(a){tapply(a, metadata_kc$treat, mean)}))

kc.PDS <- na.omit(as.data.frame(results(dds, contrast = c("treat", "PDS","con" ), cooksCutoff = FALSE))) #cooksCutoff = FALSE参数是指离群的基因不会用NA表示，之前没设置参数时，将会把离群值设为NA
kc.PDS$con <- normalized_mean_counts[match(rownames(kc.PDS),rownames(normalized_mean_counts)), 1]
kc.PDS$PDS <- normalized_mean_counts[match(rownames(kc.PDS),rownames(normalized_mean_counts)), 2]
kc.PDS$group <- ifelse(kc.PDS$pvalue<0.05&abs(kc.PDS$log2FoldChange)>=0.5,ifelse(kc.PDS$log2FoldChange>0.5,"Up","Down"),"No-sig")
table(kc.PDS$group)
# Down No-sig     Up 
# 113   8326     15
kc.PDS$geneid <- rownames(kc.PDS)
write.table(kc.PDS,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"DE.kc.PDS.txt"),
                         sep = '\t',col.names = T,row.names = F,quote = F)

kc.Phen <- na.omit(as.data.frame(results(dds, contrast = c("treat", "Phen","con"), cooksCutoff = FALSE))) #cooksCutoff = FALSE参数是指离群的基因不会用NA表示，之前没设置参数时，将会把离群值设为NA
kc.Phen$con <- normalized_mean_counts[match(rownames(kc.Phen),rownames(normalized_mean_counts)), 1]
kc.Phen$Phen <- normalized_mean_counts[match(rownames(kc.Phen),rownames(normalized_mean_counts)), 3]
kc.Phen$group <- ifelse(kc.Phen$pvalue<0.05&abs(kc.Phen$log2FoldChange)>=0.5,ifelse(kc.Phen$log2FoldChange>0.5,"Up","Down"),"No-sig")
table(kc.Phen$group)
# Down No-sig     Up 
# 224   9648    722 
kc.Phen$geneid <- rownames(kc.Phen)
write.table(kc.Phen,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"DE.kc.Phen.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)

#### Kc火山图 ####
#火山图可视化
library(ggplot2)
ggplot(kc.PDS,aes(x=log2FoldChange,y=-log10(pvalue))) +#黑白主题
  geom_point(aes(color=group),alpha=0.5,size=2) +
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
                     )+
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 0.8) +
  geom_vline(xintercept = c(-1,1),lty = 'dashed',size = 0.8) +
  annotate(geom="text", x=-1.5, y=12, size=4,label=paste0("Down\n(N=",sum(kc.PDS$group=="Down"),")"))+
  annotate(geom="text", x=1.5, y=12, size=4,label=paste0("Up\n(N=",sum(kc.PDS$group=="Up"),")"))+
  guides(color="none") + #图例颜色
  ## 修改坐标轴
  xlab(bquote(Log[2]~Fold~Change))+
  ylab(bquote(Log[10]~P~Value(PDS/DMSO))) +
  coord_cartesian(ylim = c(0,15),xlim = c(-3,2.5)) +
  labs(title="Kc167 cell") 
# ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"DEKc167vsS2.VolcanoPlot.pdf"),
       device = "pdf",width = 4,height = 3)  

ggplot(kc.Phen,aes(x=log2FoldChange,y=-log10(pvalue))) +#黑白主题
  geom_point(aes(color=group),alpha=0.5,size=2) +
  theme_bw(base_size = 16)+
  theme(#aspect.ratio = 1,
    plot.title = element_text(size = 14),
    axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12),
    panel.grid.major = element_blank(),  # 去除主要网格
    panel.grid.minor = element_blank(),  # 去除次要网格
    panel.background = element_rect(fill = "white"))+
  scale_color_manual(name = '',
                     values = c('Up'='#D6604D','Nosig'='grey','Down'='#74ADD1'), #手动设置颜色时调整颜色的因子顺序
  )+
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 0.8) +
  geom_vline(xintercept = c(-1,1),lty = 'dashed',size = 0.8) +
  annotate(geom="text", x=-3.5, y=60, size=4,label=paste0("Down\n(N=",sum(kc.Phen$group=="Down"),")"))+
  annotate(geom="text", x=4.5, y=60, size=4,label=paste0("Up\n(N=",sum(kc.Phen$group=="Up"),")"))+
  guides(color="none") + #图例颜色
  ## 修改坐标轴
  xlab(bquote(Log[2]~Fold~Change))+
  ylab(bquote(Log[10]~P~Value(PhenDC3/DMSO))) +
  coord_cartesian(ylim = c(0,120),xlim = c(-5,6)) + 
  labs(title="Kc167 cell") +
  annotate(geom = "text", x=4.7, y=115,size=4,
           label=paste0("Total:",nrow(kc.Phen))) 
# ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"DEKc167vsS2.VolcanoPlot.pdf"),
device = "pdf",width = 4,height = 3)  

gene.bed <- fread("/home/yuss/flyG4/data/ref/gene.sort.bed") %>% as.data.frame()
kc.PDS$chr <- gene.bed[match(rownames(kc.PDS),gene.bed$V4),1]
kc.PDS$start <- gene.bed[match(rownames(kc.PDS),gene.bed$V4),2]
kc.PDS$end <- gene.bed[match(rownames(kc.PDS),gene.bed$V4),3]
kc.PDS$geneid <- rownames(kc.PDS)
kc.PDS.Up <- kc.PDS[kc.PDS$group=="Up",11:13]
kc.PDS.Down <- kc.PDS[kc.PDS$group=="Down",11:13]
kc.PDS.Nosig <- kc.PDS[kc.PDS$group=="No-sig",11:13]
write.table(kc.PDS.Up,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"DE.kc.PDSUp.bed"),
sep = '\t',col.names = F,row.names = F,quote = F)
write.table(kc.PDS.Down,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"DE.kc.PDSDown.bed"),
            sep = '\t',col.names = F,row.names = F,quote = F)
write.table(kc.PDS.Nosig,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"DE.kc.PDSNosig.bed"),
            sep = '\t',col.names = F,row.names = F,quote = F)

kc.Phen$chr <- gene.bed[match(rownames(kc.Phen),gene.bed$V4),1]
kc.Phen$start <- gene.bed[match(rownames(kc.Phen),gene.bed$V4),2]
kc.Phen$end <- gene.bed[match(rownames(kc.Phen),gene.bed$V4),3]
kc.Phen.Up <- kc.Phen[kc.Phen$group=="Up",c(11:13,10)]
kc.Phen.Down <- kc.Phen[kc.Phen$group=="Down",c(11:13,10)]
kc.Phen.Nosig <- kc.Phen[kc.Phen$group=="No-sig",c(11:13,10)]
write.table(kc.Phen.Up,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"DE.kc.PhenUp.bed"),
            sep = '\t',col.names = F,row.names = F,quote = F)
write.table(kc.Phen.Down,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"DE.kc.PhenDown.bed"),
            sep = '\t',col.names = F,row.names = F,quote = F)
write.table(kc.Phen.Nosig,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"DE.kc.PhenNosig.bed"),
            sep = '\t',col.names = F,row.names = F,quote = F)

##chr X 火山图
kc.Phen.X <- kc.Phen[kc.Phen$chr=="X",]
ggplot(kc.Phen.X,aes(x=log2FoldChange,y=-log10(pvalue))) +#黑白主题
  geom_point(aes(color=group),alpha=0.5,size=2) +
  theme_bw(base_size = 16)+
  theme(#aspect.ratio = 1,
    plot.title = element_text(size = 14),
    axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12),
    panel.grid.major = element_blank(),  # 去除主要网格
    panel.grid.minor = element_blank(),  # 去除次要网格
    panel.background = element_rect(fill = "white"))+
  scale_color_manual(name = '',
                     values = c('Up'='#D6604D','Nosig'='grey','Down'='#74ADD1'), #手动设置颜色时调整颜色的因子顺序
  )+
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 0.8) +
  geom_vline(xintercept = c(-1,1),lty = 'dashed',size = 0.8) +
  annotate(geom="text", x=-1.9, y=60, size=4,label=paste0("Down\n(N=",sum(kc.Phen.X$group=="Down"),")"))+
  annotate(geom="text", x=2.7, y=60, size=4,label=paste0("Up\n(N=",sum(kc.Phen.X$group=="Up"),")"))+
  guides(color="none") + #图例颜色
  ## 修改坐标轴
  xlab(bquote(Log[2]~Fold~Change))+
  ylab(bquote(Log[10]~P~Value(PhenDC3/DMSO))) +
  coord_cartesian(ylim = c(0,120),xlim = c(-3,4)) + 
  labs(title="Kc167 cell (X chromosome)") +
  annotate(geom = "text", x=3.1, y=116,size=4,
           label=paste0("Total:",nrow(kc.Phen))) 
table(kc.Phen.X$group)

#### S2 DESeq2 ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "010.1."
readcounts_s2 <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.s2.counts.txt") %>% as.data.frame()
rownames(readcounts_s2) <- readcounts_s2$geneid
readcounts_s2 <- readcounts_s2[,-1]
metadata_s2 <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.metadata_s2.txt") %>% as.data.frame()
metadata_s2$treat <- as.factor(metadata_s2$treat)

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = readcounts_s2,
                              colData = metadata_s2,
                              design = ~ treat)
##2.dds标准化
dds <- DESeq(dds)
##3.获取标准化后的数据
normalized_counts <- counts(dds,normalized=T)
##4.有三个重复取平均
normalized_mean_counts = t(apply(normalized_counts, 1, function(a){tapply(a, metadata_s2$treat, mean)}))
# write.table(normalized_count_table,"/home/yuss/Project/Celegans_aging/RNAseq/3.Merge_result/results/normalized_counts")
# write.table(normalized_mean_count,"/home/yuss/Project/Celegans_aging/RNAseq/3.Merge_result/results/normalized_mean_counts")
# save(metadata,normalized_count_table,normalized_mean_count,dds_main1,file = "./codes/DESeq2_normalized_counts.RData")

s2.PDS <- na.omit(as.data.frame(results(dds, contrast = c("treat", "PDS", "con"), cooksCutoff = FALSE))) #cooksCutoff = FALSE参数是指离群的基因不会用NA表示，之前没设置参数时，将会把离群值设为NA
s2.PDS$con <- normalized_mean_counts[match(rownames(s2.PDS),rownames(normalized_mean_counts)), 1]
s2.PDS$PDS <- normalized_mean_counts[match(rownames(s2.PDS),rownames(normalized_mean_counts)), 2]
s2.PDS$group <- ifelse(s2.PDS$pvalue<0.05&abs(s2.PDS$log2FoldChange)>=0.5,ifelse(s2.PDS$log2FoldChange>0.5,"Up","Down"),"No-sig")
table(s2.PDS$group)
# Down No-sig     Up 
# 76   8155     44 
s2.PDS$geneid <- rownames(s2.PDS)
write.table(s2.PDS,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"DE.s2.PDS.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)

s2.Phen <- na.omit(as.data.frame(results(dds, contrast = c("treat", "Phen", "con"), cooksCutoff = FALSE))) #cooksCutoff = FALSE参数是指离群的基因不会用NA表示，之前没设置参数时，将会把离群值设为NA
s2.Phen$con <- normalized_mean_counts[match(rownames(s2.Phen),rownames(normalized_mean_counts)), 1]
s2.Phen$Phen <- normalized_mean_counts[match(rownames(s2.Phen),rownames(normalized_mean_counts)), 3]
s2.Phen$group <- ifelse(s2.Phen$pvalue<0.05&abs(s2.Phen$log2FoldChange)>=0.5,ifelse(s2.Phen$log2FoldChange>0.5,"Up","Down"),"No-sig")
table(s2.Phen$group)
# Down No-sig     Up 
# 776  10030   1351

s2.Phen$geneid <- rownames(s2.Phen)
write.table(s2.Phen,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"DE.s2.Phen.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)

#### S2火山图 ####
ggplot(s2.PDS,aes(x=log2FoldChange,y=-log10(pvalue))) +#黑白主题
  geom_point(aes(color=group),alpha=0.5,size=2) +
  theme_bw(base_size = 16)+
  theme(#aspect.ratio = 1,
    plot.title = element_text(size = 14),
    axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12),
    panel.grid.major = element_blank(),  # 去除主要网格
    panel.grid.minor = element_blank(),  # 去除次要网格
    panel.background = element_rect(fill = "white"))+
  scale_color_manual(name = '',
                     values = c('Up'='#D6604D','Nosig'='grey','Down'='#74ADD1'), #手动设置颜色时调整颜色的因子顺序
  )+
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 0.8) +
  geom_vline(xintercept = c(-1,1),lty = 'dashed',size = 0.8) +
  annotate(geom="text", x=-1.5, y=12, size=4,label=paste0("Down\n(N=",sum(s2.PDS$group=="Down"),")"))+
  annotate(geom="text", x=1.5, y=12, size=4,label=paste0("Up\n(N=",sum(s2.PDS$group=="Up"),")"))+
  guides(color="none") + #图例颜色
  ## 修改坐标轴
  xlab(bquote(Log[2]~Fold~Change))+
  ylab(bquote(Log[10]~P~Value(PDS/DMSO))) +
  coord_cartesian(ylim = c(0,20),xlim = c(-2,2)) +
  labs(title="S2 cell") 
# ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"DEs2167vsS2.VolcanoPlot.pdf"),
device = "pdf",width = 4,height = 3)  

ggplot(s2.Phen,aes(x=log2FoldChange,y=-log10(pvalue))) +#黑白主题
  geom_point(aes(color=group),alpha=0.5,size=2) +
  theme_bw(base_size = 16)+
  theme(#aspect.ratio = 1,
    plot.title = element_text(size=14), 
    axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12),
    panel.grid.major = element_blank(),  # 去除主要网格
    panel.grid.minor = element_blank(),  # 去除次要网格
    panel.background = element_rect(fill = "white"))+
  scale_color_manual(name = '',
                     values = c('Up'='#D6604D','Nosig'='grey','Down'='#74ADD1'), #手动设置颜色时调整颜色的因子顺序
  )+
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 0.8) +
  geom_vline(xintercept = c(-1,1),lty = 'dashed',size = 0.8) +
  annotate(geom="text", x=-4, y=100, size=4,label=paste0("Down\n(N=",sum(s2.Phen$group=="Down"),")"))+
  annotate(geom="text", x=6, y=100, size=4,label=paste0("Up\n(N=",sum(s2.Phen$group=="Up"),")"))+
  guides(color="none") + #图例颜色
  ## 修改坐标轴
  xlab(bquote(Log[2]~Fold~Change))+
  ylab(bquote(Log[10]~P~Value(PhenDC3/DMSO))) +
  labs(title="S2 cell") +
  coord_cartesian(ylim = c(0,200),xlim = c(-7,10)) +
  annotate(geom = "text", x=8, y=195,size=4,
           label=paste0("Total:",nrow(s2.Phen))) 
     
# ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"DEKc167vsS2.VolcanoPlot.pdf"),
device = "pdf",width = 4,height = 3)  

gene.bed <- fread("/home/yuss/flyG4/data/ref/gene.sort.bed") %>% as.data.frame()
s2.PDS$chr <- gene.bed[match(rownames(s2.PDS),gene.bed$V4),1]
s2.PDS$start <- gene.bed[match(rownames(s2.PDS),gene.bed$V4),2]
s2.PDS$end <- gene.bed[match(rownames(s2.PDS),gene.bed$V4),3]
s2.PDS.Up <- s2.PDS[s2.PDS$group=="Up",c(11:13,10)]
s2.PDS.Down <- s2.PDS[s2.PDS$group=="Down",c(11:13,10)]
s2.PDS.Nosig <- s2.PDS[s2.PDS$group=="No-sig",c(11:13,10)]
write.table(s2.PDS.Up,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"DE.s2.PDSUp.bed"),
            sep = '\t',col.names = F,row.names = F,quote = F)
write.table(s2.PDS.Down,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"DE.s2.PDSDown.bed"),
            sep = '\t',col.names = F,row.names = F,quote = F)
write.table(s2.PDS.Nosig,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"DE.s2.PDSNosig.bed"),
            sep = '\t',col.names = F,row.names = F,quote = F)

s2.Phen$chr <- gene.bed[match(rownames(s2.Phen),gene.bed$V4),1]
s2.Phen$start <- gene.bed[match(rownames(s2.Phen),gene.bed$V4),2]
s2.Phen$end <- gene.bed[match(rownames(s2.Phen),gene.bed$V4),3]
s2.Phen.Up <- s2.Phen[s2.Phen$group=="Up",c(11:13,10)]
s2.Phen.Down <- s2.Phen[s2.Phen$group=="Down",c(11:13,10)]
s2.Phen.Nosig <- s2.Phen[s2.Phen$group=="No-sig",c(11:13,10)]
write.table(s2.Phen.Up,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"DE.s2.PhenUp.bed"),
            sep = '\t',col.names = F,row.names = F,quote = F)
write.table(s2.Phen.Down,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"DE.s2.PhenDown.bed"),
            sep = '\t',col.names = F,row.names = F,quote = F)
write.table(s2.Phen.Nosig,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"DE.s2.PhenNosig.bed"),
            sep = '\t',col.names = F,row.names = F,quote = F)
##X
s2.Phen.X <- s2.Phen[s2.Phen$chr=="X",]
ggplot(s2.Phen.X,aes(x=log2FoldChange,y=-log10(pvalue))) +#黑白主题
  geom_point(aes(color=group),alpha=0.5,size=2) +
  theme_bw(base_size = 16)+
  theme(#aspect.ratio = 1,
    plot.title = element_text(size=14), 
    axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12),
    panel.grid.major = element_blank(),  # 去除主要网格
    panel.grid.minor = element_blank(),  # 去除次要网格
    panel.background = element_rect(fill = "white"))+
  scale_color_manual(name = '',
                     values = c('Up'='#D6604D','Nosig'='grey','Down'='#74ADD1'), #手动设置颜色时调整颜色的因子顺序
  )+
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 0.8) +
  geom_vline(xintercept = c(-1,1),lty = 'dashed',size = 0.8) +
  annotate(geom="text", x=-4, y=100, size=4,label=paste0("Down\n(N=",sum(s2.Phen.X$group=="Down"),")"))+
  annotate(geom="text", x=6, y=100, size=4,label=paste0("Up\n(N=",sum(s2.Phen.X$group=="Up"),")"))+
  guides(color="none") + #图例颜色
  ## 修改坐标轴
  xlab(bquote(Log[2]~Fold~Change))+
  ylab(bquote(Log[10]~P~Value(PhenDC3/DMSO))) +
  labs(title="S2 cell (X chromosome)") +
  coord_cartesian(ylim = c(0,200),xlim = c(-7,10)) +
annotate(geom = "text", x=8, y=195,size=4,
         label=paste0("Total:",nrow(s2.Phen.X))) 

s2.Phen.Up.X <- s2.Phen.Up[s2.Phen.Up$chr=="X",]
s2.Phen.Down.X <- s2.Phen.Down[s2.Phen.Down$chr=="X",]
s2.Phen.Nosig.X <- s2.Phen.Nosig[s2.Phen.Nosig$chr=="X",]
write.table(s2.Phen.Up.X,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"ChrX.DE.s2.PhenUp.bed"),
            sep = '\t',col.names = F,row.names = F,quote = F)
write.table(s2.Phen.Down.X,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"ChrX.DE.s2.PhenDown.bed"),
            sep = '\t',col.names = F,row.names = F,quote = F)
write.table(s2.Phen.Nosig.X,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"ChrX.DE.s2.PhenNosig.bed"),
            sep = '\t',col.names = F,row.names = F,quote = F)

##fisher test
data <- matrix((c(722,224,1351,776)),nrow = 2)
fisher.test(data) # 1.551e-12
data <- matrix((c(722,224,1351,776)),nrow = 2)
fisher.test(data) # 1.551e-12



#### con中Kc vs S2####
rm(list = ls());gc();rm(list = ls())#清空
Num = "010.1."
library(data.table)
readcounts <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.counts.txt") %>% as.data.frame()
metadata <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.metadata.txt") %>% as.data.frame()
metadata_con <- metadata[c(1:3,9:11),]
readcounts_con <- readcounts[,c(1:3,9:11)]
rownames(readcounts_con) <- readcounts$geneid

library(DESeq2)
##1.构建矩阵
metadata_con$cell <- as.factor(metadata_con$cell)
dds <- DESeqDataSetFromMatrix(countData = readcounts_con,
                              colData = metadata_con,
                              design = ~ cell)
##2.dds标准化
dds <- DESeq(dds)
##3.获取标准化后的数据
normalized_counts <- counts(dds,normalized=T)
##4.有三个重复取平均
normalized_mean_counts = t(apply(normalized_counts, 1, function(a){tapply(a, metadata_con$cell, mean)}))
# write.table(normalized_count_table,"/home/yuss/Project/Celegans_aging/RNAseq/3.Merge_result/results/normalized_counts")
# write.table(normalized_mean_count,"/home/yuss/Project/Celegans_aging/RNAseq/3.Merge_result/results/normalized_mean_counts")
# save(metadata,normalized_count_table,normalized_mean_count,dds_main1,file = "./codes/DESeq2_normalized_counts.RData")

con <- na.omit(as.data.frame(results(dds, contrast = c("cell", "kc","s2" ), cooksCutoff = FALSE))) #cooksCutoff = FALSE参数是指离群的基因不会用NA表示，之前没设置参数时，将会把离群值设为NA
con$kc <- normalized_mean_counts[match(rownames(con),rownames(normalized_mean_counts)), 1]
con$s2 <- normalized_mean_counts[match(rownames(con),rownames(normalized_mean_counts)), 2]
con$group <- ifelse(con$pvalue<0.05&abs(con$log2FoldChange)>=0.5,ifelse(con$log2FoldChange>0.5,"Up","Down"),"No-sig")
table(con$group)
# Down No-sig     Up 
# 3790   6298   4022  
con$geneid <- rownames(con)
write.table(con,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"DE.con.kcvss2.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)

#### phen中Kc vs S2####
rm(list = ls());gc();rm(list = ls())#清空
Num = "010.1."
readcounts <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.counts.txt") %>% as.data.frame()
metadata <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.metadata.txt") %>% as.data.frame()
metadata_treatPhen <- metadata[c(7:8,15:16),]
readcounts_treatPhen <- readcounts[,c(7:8,15:16)]
rownames(readcounts_treatPhen) <- readcounts$geneid
library(DESeq2)
##1.构建矩阵
metadata_treatPhen$cell <- as.factor(metadata_treatPhen$cell)
dds <- DESeqDataSetFromMatrix(countData = readcounts_treatPhen,
                              colData = metadata_treatPhen,
                              design = ~ cell)
##2.dds标准化
dds <- DESeq(dds)
##3.获取标准化后的数据
normalized_counts <- counts(dds,normalized=T)
##4.有三个重复取平均
normalized_mean_counts = t(apply(normalized_counts, 1, function(a){tapply(a, metadata_treatPhen$cell, mean)}))

Phen <- na.omit(as.data.frame(results(dds, contrast = c("cell", "kc","s2" ), cooksCutoff = FALSE))) #cooksCutoff = FALSE参数是指离群的基因不会用NA表示，之前没设置参数时，将会把离群值设为NA
Phen$kc <- normalized_mean_counts[match(rownames(Phen),rownames(normalized_mean_counts)), 1]
Phen$s2 <- normalized_mean_counts[match(rownames(Phen),rownames(normalized_mean_counts)), 2]
Phen$group <- ifelse(Phen$pvalue<0.05&abs(Phen$log2FoldChange)>=0.5,ifelse(Phen$log2FoldChange>0.5,"Up","Down"),"No-sig")
table(Phen$group)
# Down No-sig     Up 
# 3538   6813   3723
Phen$geneid <- rownames(Phen)
write.table(Phen,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"DE.Phen.kcvss2.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)
Phen.Up <- Phen[Phen$group=="Up",]
write.table(Phen.Up,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"DE.Phen.kcvss2Up.txt"),
            sep = '\t',col.names = F,row.names = F,quote = F)
Phen.Down <- Phen[Phen$group=="Down",]
write.table(Phen.Down,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"DE.Phen.kcvss2Down.txt"),
            sep = '\t',col.names = F,row.names = F,quote = F)

rm(list = ls());gc();rm(list = ls())#清空
Num = "010.1."
con <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.DE.con.kcvss2.txt") %>% as.data.frame()
Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.DE.Phen.kcvss2.txt") %>% as.data.frame()
con.up <- con[con$group=="Up",]
con.down <- con[con$group=="Down",]
Phen.up <- Phen[Phen$group=="Up",]
Phen.down <- Phen[Phen$group=="Down",]
upresult <- left_join(con.up, Phen.up, by = "geneid")
upresult <- na.omit(left_join(con.up, Phen.up, by = "geneid")) #3216
upgene <- upresult[upresult$log2FoldChange.x<upresult$log2FoldChange.y, ] #1632
uplist <- unlist(upgene$geneid)

##富集
library(clusterProfiler)
library(HDO.db)
library(DOSE)
library(org.Dm.eg.db)
enrich_go_kegg<-function(gene_id){
  result<-list()
  ensembl_2_entrezid<-bitr(gene_id,fromType ="ENSEMBL", toType = c("ENTREZID","SYMBOL"), OrgDb = org.Dm.eg.db)
  ego_ALL <- enrichGO(gene = as.character(ensembl_2_entrezid$ENTREZID),
                      OrgDb=org.Dm.eg.db,
                      keyType = "ENTREZID",
                      ont = "ALL",
                      pAdjustMethod = "BH",
                      minGSSize = 1,
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05,
                      readable = TRUE)
  result$ego_all<-as.data.frame(ego_ALL)
  ensembl_2_kegg_id<-bitr_kegg(ensembl_2_entrezid$ENTREZID,fromType = "ncbi-geneid",toType = "kegg",organism = "dme")
  kegg<-enrichKEGG(gene = ensembl_2_entrezid$ENTREZID,organism = "dme", keyType = "ncbi-geneid",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
  result$kegg<-kegg@result
  result$ego_ALL<-ego_ALL
  return(result)
}

enrich_up <- enrich_go_kegg(uplist)
enrich_up_ego <- enrich_up$ego_all
enrich_up_kegg <- enrich_up$kegg
enrich_up_ego_ALL <- enrich_up$ego_ALL
dotplot(enrich_up_ego_ALL,showCategory=10,split="ONTOLOGY",label_format=60,title="Up regulated GO enrichment") + facet_grid(ONTOLOGY~., scale='free')
##showCategory指定展示的GO Terms的个数，默认展示显著富集的top10个
##label_format=60左边的名称一行显示60
##facet_grid(ONTOLOGY~., scale='free')将三个基因本体分开
library(ggplot2)
enrich_up_kegg <- arrange(enrich_up_kegg, enrich_up_kegg$p.adjust)
enrich_up_kegg$Description1 <- sapply(strsplit(enrich_up_kegg$Description, " - "), function(x) x[1])
ggplot(enrich_up_kegg[1:14,],aes(x=Count/116,y=Description1,colour=-1*log10(pvalue),size=Count))+
  geom_point()+
  scale_size(range = c(2, 8))+ #scale_size修改图中的点的大小，范围是2到8
  scale_color_gradient(low = "blue", high = "red")+
  theme_bw()+
  ylab("KEGG Pathway Terms")+
  xlab("Gene Ratio")+
  labs(color=expression(-log[10](PValue)), title = "Up regulated kegg enrichment") +#expression函数改变样式，[]是用来添加下标，^是用来添加上标
  theme(text = element_text(size=15))


upresult <- na.omit(left_join(con.up, Phen.down, by = "geneid")) #68
downresult <- na.omit(left_join(con.down, Phen.up, by = "geneid")) #74
downresult <- na.omit(left_join(con.down, Phen.down, by = "geneid")) #3099
downgene <- downresult[downresult$log2FoldChange.x>downresult$log2FoldChange.y, ] #1176
downlist <- unlist(downgene$geneid)

enrich_down <- enrich_go_kegg(downlist)
enrich_down_ego <- enrich_down$ego_all
enrich_down_kegg <- enrich_down$kegg
enrich_down_ego_ALL <- enrich_down$ego_ALL
dotplot(enrich_down_ego_ALL,showCategory=10,split="ONTOLOGY",label_format=60,title="Down regulated GO enrichment") + facet_grid(ONTOLOGY~., scale='free')
##showCategory指定展示的GO Terms的个数，默认展示显著富集的top10个
##label_format=60左边的名称一行显示60
##facet_grid(ONTOLOGY~., scale='free')将三个基因本体分开
library(ggplot2)
enrich_up_kegg <- arrange(enrich_up_kegg, enrich_up_kegg$p.adjust)
enrich_up_kegg$Description1 <- sapply(strsplit(enrich_up_kegg$Description, " - "), function(x) x[1])
ggplot(enrich_up_kegg[1:14,],aes(x=Count/116,y=Description1,colour=-1*log10(pvalue),size=Count))+
  geom_point()+
  scale_size(range = c(2, 8))+ #scale_size修改图中的点的大小，范围是2到8
  scale_color_gradient(low = "blue", high = "red")+
  theme_bw()+
  ylab("KEGG Pathway Terms")+
  xlab("Gene Ratio")+
  labs(color=expression(-log[10](PValue)), title = "Up regulated kegg enrichment") +#expression函数改变样式，[]是用来添加下标，^是用来添加上标
  theme(text = element_text(size=15))