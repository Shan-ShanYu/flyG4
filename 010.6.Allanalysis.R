#*差异表达分析---------------------------------------------------------------------
#### 1.Kc Phen vs con DESeq2 #### 
rm(list = ls());gc();rm(list = ls())#清空
Num = "010.6."
readcounts_kc <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.kc.counts.txt") %>% as.data.frame()
rownames(readcounts_kc) <- readcounts_kc$geneid
readcounts_kc <- readcounts_kc[,-c(1,5:7)]
metadata_kc <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.metadata_kc.txt") %>% as.data.frame()
metadata_kc <- metadata_kc[-c(4:6),]
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

kc.Phen <- na.omit(as.data.frame(results(dds, contrast = c("treat", "Phen","con"), cooksCutoff = FALSE))) #cooksCutoff = FALSE参数是指离群的基因不会用NA表示，之前没设置参数时，将会把离群值设为NA
kc.Phen$con <- normalized_mean_counts[match(rownames(kc.Phen),rownames(normalized_mean_counts)), 1]
kc.Phen$Phen <- normalized_mean_counts[match(rownames(kc.Phen),rownames(normalized_mean_counts)), 2]
kc.Phen$group <- ifelse(kc.Phen$pvalue<0.05&abs(kc.Phen$log2FoldChange)>=0.5,ifelse(kc.Phen$log2FoldChange>0.5,"Up","Down"),"No-sig")
table(kc.Phen$group)
# Down No-sig     Up 
# 211   8967    655 
kc.Phen$geneid <- rownames(kc.Phen)
write.table(kc.Phen,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"DE.kc.Phen.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)

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
  geom_vline(xintercept = c(-0.5,0.5),lty = 'dashed',size = 0.8) +
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
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"DEKc.Phenvscon.VolcanoPlot.pdf"),
device = "pdf",width = 4,height = 3)  


#### 2.S2 Phen vs con DESeq2 #### 
readcounts_s2 <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.s2.counts.txt") %>% as.data.frame()
rownames(readcounts_s2) <- readcounts_s2$geneid
readcounts_s2 <- readcounts_s2[,-c(1,5:7)]
metadata_s2 <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.metadata_s2.txt") %>% as.data.frame()
metadata_s2 <- metadata_s2[-c(4:6),]
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

s2.Phen <- na.omit(as.data.frame(results(dds, contrast = c("treat", "Phen", "con"), cooksCutoff = FALSE))) #cooksCutoff = FALSE参数是指离群的基因不会用NA表示，之前没设置参数时，将会把离群值设为NA
s2.Phen$con <- normalized_mean_counts[match(rownames(s2.Phen),rownames(normalized_mean_counts)), 1]
s2.Phen$Phen <- normalized_mean_counts[match(rownames(s2.Phen),rownames(normalized_mean_counts)), 2]
s2.Phen$group <- ifelse(s2.Phen$pvalue<0.05&abs(s2.Phen$log2FoldChange)>=0.5,ifelse(s2.Phen$log2FoldChange>0.5,"Up","Down"),"No-sig")
table(s2.Phen$group)
# Down No-sig     Up 
# 740   9865   1306
s2.Phen$geneid <- rownames(s2.Phen)
write.table(s2.Phen,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"DE.s2.Phen.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)

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
                     values = c('Up'='#BF69A9','Nosig'='grey','Down'='#8BC25F'), #手动设置颜色时调整颜色的因子顺序
  )+
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 0.8) +
  geom_vline(xintercept = c(-0.5,0.5),lty = 'dashed',size = 0.8) +
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
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"DES2.Phenvscon.VolcanoPlot.pdf"),
       device = "pdf",width = 4,height = 3)  

#### 3.共同的上下调基因 ####
# 筛选上调的基因
upregulated_genes <- intersect(kc.Phen$geneid[kc.Phen$group=="Up"], s2.Phen$geneid[s2.Phen$group=="Up"])

# 筛选下调的基因
downregulated_genes <- intersect(kc.Phen$geneid[kc.Phen$group=="Down"], s2.Phen$geneid[s2.Phen$group=="Down"])

# 统计数量
length(upregulated_genes) #共同上调299
length(downregulated_genes) #34

library(eulerr)
VennDiag <- euler(c("kcup" = 356,"s2up" = 1007,
                    "kcup&s2up" = 299))
p <- plot(VennDiag, counts = FALSE, font=3, cex=1, alpha=1,quantities = TRUE,lwd =3.5,
     labels=c("Kc167 Up","S2 Up"),
     label.col = "white",fill="white",
     col = c('#D6604D','#BF69A9'))
ggsave(p,filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"KcS2UpDEG.Intersect.pdf"),
       device = "pdf",width = 3.8,height = 3)
library(eulerr)
VennDiag <- euler(c("kcdown" = 177,"s2down" = 706,
                    "kcdown&s2down" = 34))
p <- plot(VennDiag, counts = FALSE, font=3, cex=1, alpha=1,quantities = TRUE,lwd =3.5,
     labels=c("Kc167 Down","S2 Down"),
     label.col = "white",fill="white",
     col = c('#74ADD1','#8BC25F'))
ggsave(p,filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"KcS2DownDEG.Intersect.pdf"),
       device = "pdf",width = 3.8,height = 3)

#统计上调共有基因且含有eG4的数量
gene.kc <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.gene.kc.bed") %>% as.data.frame()
kc.Phen$num <- gene.kc[match(kc.Phen$geneid,gene.kc$V4),7]
kc.Phen$kc.eG4 <- ifelse(kc.Phen$num==0,"no eG4","eG4")
column_vector <- data.frame(upregulated_genes)
column_vector$kc.eG4 <- kc.Phen[match(column_vector$upregulated_genes,kc.Phen$geneid),12]
gene.s2 <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.gene.s2.bed") %>% as.data.frame()
s2.Phen$num <- gene.s2[match(s2.Phen$geneid,gene.s2$V4),7]
s2.Phen$s2.eG4 <- ifelse(s2.Phen$num==0,"no eG4","eG4")
column_vector$s2.eG4 <- s2.Phen[match(column_vector$upregulated_genes,s2.Phen$geneid),12]
eG4up <- column_vector[c(column_vector$s2.eG4 == "eG4" & column_vector$kc.eG4 == "eG4"), ]

#统计下调共有基因且含有eG4的数量
down_vector <- data.frame(downregulated_genes)
down_vector$kc.eG4 <- kc.Phen[match(down_vector$downregulated_genes,kc.Phen$geneid),12]
down_vector$s2.eG4 <- s2.Phen[match(down_vector$downregulated_genes,s2.Phen$geneid),12]
eG4down <- down_vector[(down_vector$s2.eG4 == "eG4" | down_vector$kc.eG4 == "eG4"), ]

library(clusterProfiler)
library(HDO.db)
library(DOSE)
library(dplyr)
library(data.table)

commonup <- unlist(eG4up$upregulated_genes)

##富集分析
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
enrich_commonup <- enrich_go_kegg(commonup)
enrich_commonup_ego <- enrich_commonup$ego_all  ##没有GO
enrich_commonup_kegg <- enrich_commonup$kegg     ##4个
enrich_commonup_ego_ALL <- enrich_commonup$ego_ALL


#上下调方向不同
upregulated_genes <- intersect(kc.Phen$geneid[kc.Phen$group=="Up"], s2.Phen$geneid[s2.Phen$group=="Down"])
downregulated_genes <- intersect(kc.Phen$geneid[kc.Phen$group=="Down"], s2.Phen$geneid[s2.Phen$group=="Up"])
# 统计数量
length(upregulated_genes) #35
length(downregulated_genes) #7

VennDiag <- euler(c("kcup" = 655-35,"s2down" = 740-35,
                    "kcup&s2down" = 35))
p <- plot(VennDiag, counts = FALSE, font=3, cex=1, alpha=1,quantities = TRUE,lwd =3.5,
          labels=c("Kc167 Up","S2 Down"),
          label.col = "white",fill="white",
          col = c('#D6604D','#8BC25F'))
p
ggsave(p,filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"KcUpS2DownDEG.Intersect.pdf"),
       device = "pdf",width = 3.8,height = 3)

VennDiag <- euler(c("kcdown" = 211-7,"s2up" = 1306-7,
                    "kcdown&s2up" = 7))
p <- plot(VennDiag, counts = FALSE, font=3, cex=1, alpha=1,quantities = TRUE,lwd =3.5,
          labels=c("Kc167 Down","S2 Up"),
          label.col = "white",fill="white",
          col = c('#74ADD1','#BF69A9'))
p
ggsave(p,filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"KcDownS2UpDEG.Intersect.pdf"),
       device = "pdf",width = 3.8,height = 3)

#*eG4信号------------------------------------------------------------------------------
#### 1.Kc fisher #### 
rm(list = ls());gc();rm(list = ls())#清空
Num = "010.6."
kc.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.6.DE.kc.Phen.txt") %>% as.data.frame()
gene.kc <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.gene.kc.bed") %>% as.data.frame()
kc.Phen$num <- gene.kc[match(kc.Phen$geneid,gene.kc$V4),7]
kc.Phen$kc.eG4 <- ifelse(kc.Phen$num==0,"no eG4","eG4")

table(kc.Phen$kc.eG4,kc.Phen$group)
# Down No-sig   Up
# eG4      29   2100  300
# no eG4  182   6867  355
# DEG eG4=29+300=329
# DEG no eG4=182+355=537
data <- matrix(round((c(329,2100,537,6867)),0),nrow = 2)
fisher.test(data) # p-value < 2.2e-16
data <- matrix(round((c(300,2100,355,6867)),0),nrow = 2)
fisher.test(data) # p-value < 2.2e-16
data <- matrix(round((c(29,2100,182,6867)),0),nrow = 2) #下调基因缺失eG4
fisher.test(data) # p-value < 2.2e-16

#### Kc基因上eG4的信号####
kc.Phen$chr <- gene.kc[match(kc.Phen$geneid,gene.kc$V4),1]
kc.Phen$start <- gene.kc[match(kc.Phen$geneid,gene.kc$V4),2]
kc.Phen$end <- gene.kc[match(kc.Phen$geneid,gene.kc$V4),3]
kc.Phen$length <- kc.Phen$end-kc.Phen$start
kc.Phen$density <- kc.Phen$num/kc.Phen$length
my_comparisons = list(c("No-sig","Down"),c("No-sig","Up"))
ggplot(data = kc.Phen,aes(x=group,y=density*1000,fill=group)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(0.45,0.55),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(name = '',
                    values = c('Up'='#D6604D','No-sig'='grey','Down'='#74ADD1'), #手动设置颜色时调整颜色的因子顺序
  )+
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Number of eG4 (per kb)") +
  coord_cartesian(ylim = c(0,1)) +
  labs(title="Kc167 cell")
tapply(kc.Phen$density, kc.Phen$group, mean)
# Down       No-sig           Up 
# 4.610469e-05 7.377292e-05 1.950514e-04 
tapply(kc.Phen$density, kc.Phen$group, median)

kc.Phen$numtype <- ifelse(kc.Phen$num > 0,"have","no")
a <- as.data.frame(table(kc.Phen$group, kc.Phen$numtype))
wide_a <- spread(a, Var2, Freq)
wide_a$sum <- wide_a$have + wide_a$no
wide_a$have.ratio <- wide_a$have/wide_a$sum
wide_a$no.ratio <- wide_a$no/wide_a$sum
wide_a$Var1 <- factor(wide_a$Var1,levels = c("Down","No-sig","Up"))
  
library(ggsci)
ggplot(wide_a, aes(x=Var1,y=have.ratio*100,fill=Var1)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(0.7),width = 0.7,color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值
  coord_cartesian(ylim = c(0,60)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(name = '',
                    values = c('Up'='#D6604D','No-sig'='grey','Down'='#74ADD1'), #手动设置颜色时调整颜色的因子顺序
  )+
  # geom_text(aes(label=have.ratio*100,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题(左下边框，没有网格线)
  ylab("(%) Gene with eG4") + ##通过bquote函数给图标签添加上下标
  geom_signif(y_position=26.5, xmin=1, xmax=2,
              annotation=c("6.8e-04"),tip_length=0)+
  geom_signif(y_position=49, xmin=2, xmax=3,
              annotation=c("2.2e-16"),tip_length=0)+
  geom_signif(y_position=56, xmin=1, xmax=3,
              annotation=c("2.2e-16"),tip_length=0)+
  theme(plot.title = element_text(size = 14, face = "plain"),
        axis.title.y = element_text(size = 14), ##y坐标轴标题字体大小
        axis.title.x = element_blank(), ##删除x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.position = "none") +
  labs(title="Kc167 cell")
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"KcPhenDEG.PercentagewitheG4.pdf"),
       device = "pdf",width = 3.2,height = 3.4)
data <- matrix((c(29,2100,182,6867)),nrow = 2)
fisher.test(data) #0.0006784

data <- matrix((c(300,2100,355,6867)),nrow = 2)
fisher.test(data) #p-value < 2.2e-16

data <- matrix((c(300,29,355,182)),nrow = 2)
fisher.test(data) #p-value < 2.2e-16

#### Kc基因的Promoter区eG4信号 ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "010.6."
promoter.kcG4 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter2000.kc.bed") %>% as.data.frame()
kc.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.6.DE.kc.Phen.txt") %>% as.data.frame()
table(kc.Phen$group)
kc.Phen$promoter.num <- promoter.kcG4[match(kc.Phen$geneid,promoter.kcG4$V4),6]
kc.Phen$group <- factor(kc.Phen$group,levels = c("No-sig","Down","Up"))
my_comparisons = list(c("No-sig","Down"),c("No-sig","Up"))
ggplot(data = kc.Phen,aes(x=group,y=promoter.num,fill=group)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(3,3.5,4,4.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  # scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Number of Promoter eG4") +
  coord_cartesian(ylim = c(0,5)) +
  scale_y_continuous(expand = c(0,0)) 
tapply(kc.Phen$promoter.num, kc.Phen$group, mean)
# No-sig      Down        Up 
# 0.3038920 0.1943128 0.5221374

kc.Phen$promter.numtype <- ifelse(kc.Phen$promoter.num > 0,"have","no")
a <- as.data.frame(table(kc.Phen$group, kc.Phen$promter.numtype))
wide_a <- spread(a, Var2, Freq)
wide_a$sum <- wide_a$have + wide_a$no
wide_a$have.ratio <- wide_a$have/wide_a$sum
wide_a$no.ratio <- wide_a$no/wide_a$sum
wide_a$Var1 <- factor(wide_a$Var1,levels = c("Down","No-sig","Up"))
ggplot(wide_a, aes(x=Var1,y=have.ratio*100,fill=Var1)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(0.7),width = 0.7,color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值
  coord_cartesian(ylim = c(0,40)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(name = '',
                    values = c('Up'='#D6604D','No-sig'='grey','Down'='#74ADD1'), #手动设置颜色时调整颜色的因子顺序
  )+
  # geom_text(aes(label=have.ratio*100,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题(左下边框，没有网格线)
  ylab("(%) Genes with eG4 in the promoter") + ##通过bquote函数给图标签添加上下标
  geom_signif(y_position=26, xmin=1, xmax=2,
              annotation=c("0.032"),tip_length=0)+
  geom_signif(y_position=35, xmin=2, xmax=3,
              annotation=c("6.69e-10"),tip_length=0)+
  geom_signif(y_position=37.5, xmin=1, xmax=3,
              annotation=c("7.73e-07"),tip_length=0)+
  theme(plot.title = element_text(size = 14, face = "plain"),
        axis.title.y = element_text(size = 14), ##y坐标轴标题字体大小
        axis.title.x = element_blank(), ##删除x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.position = "none") +
  labs(title="Kc167 cell")
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"KcPhenDEG.PercentagewithPromotereG4.pdf"),
       device = "pdf",width = 3.2,height = 3.4)

data <- matrix((c(2107,36,6860,175)),nrow = 2)
fisher.test(data) #0.03194

data <- matrix((c(2107,227,6860,428)),nrow = 2)
fisher.test(data) #6.693e-10

data <- matrix((c(227,36,428,175)),nrow = 2)
fisher.test(data) #7.732e-07


#### S2 fisher####
rm(list = ls());gc();rm(list = ls())#清空
Num = "010.6."
s2.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.6.DE.s2.Phen.txt") %>% as.data.frame()
gene.s2 <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.gene.s2.bed") %>% as.data.frame()
s2.Phen$num <- gene.s2[match(s2.Phen$geneid,gene.s2$V4),7]
s2.Phen$s2.eG4 <- ifelse(s2.Phen$num==0,"no eG4","eG4")
table(s2.Phen$s2.eG4,s2.Phen$group)
# Down No-sig   Up
# eG4     253   2088  369
# no eG4  487   7777  937
# DEG eG4= 253+369=622
# DEG no eG4=487+937=1424
data <- matrix(round((c(622,2088,1424,7777)),0),nrow = 2)
fisher.test(data) # p-value < 2.2e-16

#### 2.S2基因上eG4的信号(密度) ####
s2.Phen$chr <- gene.s2[match(s2.Phen$geneid,gene.s2$V4),1]
s2.Phen$start <- gene.s2[match(s2.Phen$geneid,gene.s2$V4),2]
s2.Phen$end <- gene.s2[match(s2.Phen$geneid,gene.s2$V4),3]
s2.Phen$length <- s2.Phen$end-s2.Phen$start
s2.Phen$density <- s2.Phen$num/s2.Phen$length

s2.Phen$group <- factor(s2.Phen$group,levels = c("Down","No-sig","Up"))
my_comparisons = list(c("No-sig","Down"),c("No-sig","Up"),c("Down","Up"))

ggplot(data = s2.Phen,aes(x=group,y=density*1000,fill=group)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(-0.18,-0.15,-0.175),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(name = '',
                    values = c('Up'='#D6604D','No-sig'='grey','Down'='#74ADD1'), #手动设置颜色时调整颜色的因子顺序
  )+
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Number of eG4 (per kb)") +
  coord_cartesian(ylim = c(0,0.3)) +
  labs(title="S2 cell")

tapply(s2.Phen$density, s2.Phen$group, mean)
# Down       No-sig           Up 
# 9.318734e-05 7.480327e-05 1.434526e-04

s2.Phen$numtype <- ifelse(s2.Phen$num > 0,"have","no")
a <- as.data.frame(table(s2.Phen$group, s2.Phen$numtype))
wide_a <- spread(a, Var2, Freq)
wide_a$sum <- wide_a$have + wide_a$no
wide_a$have.ratio <- wide_a$have/wide_a$sum
wide_a$no.ratio <- wide_a$no/wide_a$sum

ggplot(wide_a, aes(x=Var1,y=have.ratio*100,fill=Var1)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(0.7),width = 0.7,color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值
  coord_cartesian(ylim = c(0,41)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_npg()+ ##更改颜色
  # geom_text(aes(label=have.ratio*100,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题(左下边框，没有网格线)
  ylab("(%) Gene with eG4") + ##通过bquote函数给图标签添加上下标
  # geom_signif(y_position=35, xmin=1, xmax=2,
  #             annotation=c("6.0e-14"),tip_length=0)+ 
  # geom_signif(y_position=37, xmin=2, xmax=3,
  #             annotation=c("0.0076"),tip_length=0)+ 
  # geom_signif(y_position=39, xmin=1, xmax=3,
  #             annotation=c("5.5e-08"),tip_length=0)+ 
  theme(plot.title = element_text(size = 14),
        axis.title.y = element_text(size = 16), ##y坐标轴标题字体大小
        axis.title.x = element_blank(), ##删除x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.title = element_blank()) +
  labs(title="S2 cell")

#### S2基因的Promoter区eG4的信号 ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "010.6."
promoter.s2G4 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter2000.s2.bed") %>% as.data.frame()
s2.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.6.DE.s2.Phen.txt") %>% as.data.frame()
table(s2.Phen$group)
s2.Phen$promoter.num <- promoter.s2G4[match(s2.Phen$geneid,promoter.s2G4$V4),6]
s2.Phen$group <- factor(s2.Phen$group,levels = c("Down","No-sig","Up"))
my_comparisons = list(c("No-sig","Down"),c("No-sig","Up"))
ggplot(data = s2.Phen,aes(x=group,y=promoter.num,fill=group)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(3,3.5,4,4.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  # scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Number of Promoter eG4") +
  coord_cartesian(ylim = c(0,10)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title="S2 cell")
tapply(s2.Phen$promoter.num, s2.Phen$group, mean)
# Down    No-sig        Up 
# 0.3297297 0.3049164 0.4127106

s2.Phen$promter.numtype <- ifelse(s2.Phen$promoter.num > 0,"have","no")
a <- as.data.frame(table(s2.Phen$group, s2.Phen$promter.numtype))
wide_a <- spread(a, Var2, Freq)
wide_a$sum <- wide_a$have + wide_a$no
wide_a$have.ratio <- wide_a$have/wide_a$sum
wide_a$no.ratio <- wide_a$no/wide_a$sum


ggplot(wide_a, aes(x=Var1,y=have.ratio*100,fill=Var1)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(0.7),width = 0.7,color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值
  coord_cartesian(ylim = c(0,40)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(name = '',
                    values = c('Up'='#BF69A9','No-sig'='grey','Down'='#8BC25F'), #手动设置颜色时调整颜色的因子顺序
  )+
  cowplot::theme_half_open() + ##主题(左下边框，没有网格线)
  ylab("(%) Gene with eG4 in the promoter") + ##通过bquote函数给图标签添加上下标
  geom_signif(y_position=25, xmin=1, xmax=2,
              annotation=c("0.649"),tip_length=0)+
  geom_signif(y_position=35, xmin=1, xmax=3,
              annotation=c("0.008"),tip_length=0)+
  geom_signif(y_position=31, xmin=2, xmax=3,
              annotation=c("1.73e-06"),tip_length=0)+
  theme(plot.title = element_text(size = 14, face = "plain"),
        axis.title.y = element_text(size = 14), ##y坐标轴标题字体大小
        axis.title.x = element_blank(), ##删除x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.position = "none") +
  labs(title="S2 cell")
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"s2PhenDEG.PercentagewithPromotereG4.pdf"),
       device = "pdf",width = 3.2,height = 3.4)
wide_a
data <- matrix((c(172,2224,568,7641)),nrow = 2)
fisher.test(data) # p-value =0.6489
data <- matrix((c(374,2224,932,7641)),nrow = 2)
fisher.test(data) # p-value =1.732e-06
data <- matrix((c(172,374,568,932)),nrow = 2)
fisher.test(data) # p-value =0.008012

#*差异倍数 log2FC---------------------------------------------------------------------------
#### 1.Kc 基因含有eG4和no eG4差异倍数  ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "010.6."
kc.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.6.DE.kc.Phen.txt") %>% as.data.frame()
gene.kc <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.gene.kc.bed") %>% as.data.frame()
kc.Phen$num <- gene.kc[match(kc.Phen$geneid,gene.kc$V4),7]
kc.Phen$kc.eG4 <- ifelse(kc.Phen$num==0,"no eG4","eG4")
kc.Phen$kc.eG4 <- factor(kc.Phen$kc.eG4,levels = c("no eG4","eG4"))

my_comparisons = list(c("eG4","no eG4"))
ggplot(data = kc.Phen,aes(x=kc.eG4,y=log2FoldChange,fill=kc.eG4,color=kc.eG4)) +
  geom_violin(position = position_dodge(width = 1), scale = 'width') +
  geom_boxplot(position = position_dodge(width = 1), outlier.size = 0.6, width = 0.2,
               show.legend = FALSE) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(5,5),
 tip.length = 0,size=4)+
  coord_cartesian(ylim = c(-3,6)) +
  scale_fill_manual(values = c("#818181","#D08725")) +
  scale_color_manual(values = c("#5b5b5b","#8f5d19")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2](PhenDC3/DMSO))) +
  labs(title="Kc167 cell") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_discrete(labels = c("without eG4","eG4"))
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"KcPhenlogFC.witheG4noeG4.pdf"),
       device = "pdf",width = 2.7,height = 2.7)

tapply(kc.Phen$log2FoldChange,kc.Phen$kc.eG4, mean)
tapply(kc.Phen$log2FoldChange,kc.Phen$kc.eG4, median)

sample_data <- kc.Phen[kc.Phen$kc.eG4=="eG4",2]
# 单样本Wilcoxon检验(看大于0是否有显著性)
wilcox.test(sample_data, mu = 0,alternative = "greater")
sample_data <- kc.Phen[kc.Phen$kc.eG4=="no eG4",2]
# 单样本Wilcoxon检验(看小于0是否有显著性)
wilcox.test(sample_data, mu = 0,alternative = "less") #p-value < 2.2e-16

#### 2.Kc 启动子含有eG4和no eG4差异倍数 ####
rm(list = ls());gc();rm(list = ls())
Num = "010.6."
promoter.kcG4 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter2000.kc.bed") %>% as.data.frame()
kc.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.6.DE.kc.Phen.txt") %>% as.data.frame()
kc.Phen$pro.num <- promoter.kcG4[match(kc.Phen$geneid,promoter.kcG4$V4),6]
kc.Phen$kc.eG4 <- ifelse(kc.Phen$pro.num==0,"no eG4","eG4")
kc.Phen$kc.eG4 <- factor(kc.Phen$kc.eG4,levels = c("no eG4","eG4"))

my_comparisons = list(c("eG4","no eG4"))
ggplot(data = kc.Phen,aes(x=kc.eG4,y=log2FoldChange,fill=kc.eG4,color=kc.eG4)) +
  geom_violin(position = position_dodge(width = 1), scale = 'width') +
  geom_boxplot(position = position_dodge(width = 1), outlier.size = 0.6, width = 0.2, show.legend = FALSE) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(5,5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(-3,6)) +
  scale_fill_manual(values = c("#818181","#D08725")) +
  scale_color_manual(values = c("#5b5b5b","#8f5d19")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2](PhenDC3/DMSO))) +
  labs(title="Kc167 cell (Promoter eG4)") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_discrete(labels = c("without eG4","eG4"))
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"KcPhenlogFC.withPromotereG4noeG4.pdf"),
       device = "pdf",width = 2.7,height = 2.7)

tapply(kc.Phen$log2FoldChange,kc.Phen$kc.eG4, mean)
tapply(kc.Phen$log2FoldChange,kc.Phen$kc.eG4, median)

# 单样本Wilcoxon检验(看大于0是否有显著性)
sample_data <- kc.Phen[kc.Phen$kc.eG4=="eG4",2]
wilcox.test(sample_data, mu = 0,alternative = "greater") #p-value = 2.7e-16
# 单样本Wilcoxon检验(看小于0是否有显著性)
sample_data <- kc.Phen[kc.Phen$kc.eG4=="no eG4",2]
wilcox.test(sample_data, mu = 0,alternative = "less") #p-value = 0.001124

#### 3.S2 基因含有eG4和no eG4差异倍数  ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "010.6."
s2.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.6.DE.s2.Phen.txt") %>% as.data.frame()
gene.s2 <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.gene.s2.bed") %>% as.data.frame()
s2.Phen$num <- gene.s2[match(s2.Phen$geneid,gene.s2$V4),7]
s2.Phen$s2.eG4 <- ifelse(s2.Phen$num==0,"no eG4","eG4")
s2.Phen$s2.eG4 <- factor(s2.Phen$s2.eG4,levels = c("no eG4","eG4"))
my_comparisons = list(c("eG4","no eG4"))
ggplot(data = s2.Phen,aes(x=s2.eG4,y=log2FoldChange,fill=s2.eG4,color=s2.eG4)) +
  geom_violin(position = position_dodge(width = 1), scale = 'width') +
  geom_boxplot(position = position_dodge(width = 1), outlier.size = 0.6, width = 0.2, show.legend = FALSE) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(8.5,5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(-5,10)) +
  scale_fill_manual(values = c("#818181","#D08725")) +
  scale_color_manual(values = c("#5b5b5b","#8f5d19")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2](PhenDC3/DMSO))) +
  labs(title="S2 cell") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_discrete(labels = c("without eG4","eG4"))
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"S2PhenlogFC.witheG4noeG4.pdf"),
       device = "pdf",width = 2.7,height = 2.7)

tapply(s2.Phen$log2FoldChange,s2.Phen$s2.eG4, mean)
tapply(s2.Phen$log2FoldChange,s2.Phen$s2.eG4, median)

sample_data <- s2.Phen[s2.Phen$s2.eG4=="eG4",2]
# 单样本Wilcoxon检验(看大于0是否有显著性)
wilcox.test(sample_data, mu = 0,alternative = "greater") #p-value = 1.049e-13
sample_data <- s2.Phen[s2.Phen$s2.eG4=="no eG4",2]
# 单样本Wilcoxon检验(看小于0是否有显著性)
wilcox.test(sample_data, mu = 0,alternative = "greater") #p-value = 1.423e-12 显著大于0

#### 4.S2 启动子含有eG4和no eG4差异倍数  ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "010.6."
s2.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.6.DE.s2.Phen.txt") %>% as.data.frame()
promoter.s2G4 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter2000.s2.bed") %>% as.data.frame()
s2.Phen$pro.num <- promoter.s2G4[match(s2.Phen$geneid,promoter.s2G4$V4),6]
s2.Phen$s2.eG4 <- ifelse(s2.Phen$pro.num==0,"no eG4","eG4")
s2.Phen$s2.eG4 <- factor(s2.Phen$s2.eG4,levels = c("no eG4","eG4"))
my_comparisons = list(c("eG4","no eG4"))
ggplot(data = s2.Phen,aes(x=s2.eG4,y=log2FoldChange,fill=s2.eG4,color=s2.eG4)) +
  geom_violin(position = position_dodge(width = 1), scale = 'width') +
  geom_boxplot(position = position_dodge(width = 1), outlier.size = 0.6, width = 0.2, show.legend = FALSE) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(8.5,5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(-5,10)) +
  scale_fill_manual(values = c("#818181","#D08725")) +
  scale_color_manual(values = c("#5b5b5b","#8f5d19")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2](PhenDC3/DMSO))) +
  labs(title="S2 cell (Promter eG4)") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_discrete(labels = c("without eG4","eG4"))
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"S2PhenlogFC.withPromotereG4noeG4.pdf"),
       device = "pdf",width = 2.7,height = 2.7)

tapply(s2.Phen$log2FoldChange,s2.Phen$s2.eG4, mean)
tapply(s2.Phen$log2FoldChange,s2.Phen$s2.eG4, median)

sample_data <- s2.Phen[s2.Phen$s2.eG4=="eG4",2]
# 单样本Wilcoxon检验(看大于0是否有显著性)
wilcox.test(sample_data, mu = 0,alternative = "greater") #p-value < 2.2e-16 显著大于0
sample_data <- s2.Phen[s2.Phen$s2.eG4=="no eG4",2]
# 单样本Wilcoxon检验(看小于0是否有显著性)
wilcox.test(sample_data, mu = 0,alternative = "greater") #p-value = 6.303e-09 显著大于0

#*分成三类看差异倍数------------------------------------------------------------------
#### Kc基因含有eG4\non-eG4\withouteG4差异倍数 ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "010.6."
kc.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.6.DE.kc.Phen.txt") %>% as.data.frame()
gene.pqs <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.2.gene.pqs.bed") %>% as.data.frame()
##一样的gene.pqs <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.gene.pqs.bed") %>% as.data.frame()
kc.Phen$pqsnum <- gene.pqs[match(kc.Phen$geneid,gene.pqs$V4),7]
kc.Phen$pqs.type <- ifelse(kc.Phen$pqsnum==0,"other","pqs")
gene.kc <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.gene.kc.bed") %>% as.data.frame()
kc.Phen$num <- gene.kc[match(kc.Phen$geneid,gene.kc$V4),7]
kc.Phen$kc.eG4 <- ifelse(kc.Phen$pqs.type=="pqs",ifelse(kc.Phen$num=="0","non-eG4","eG4"),"other") #增加kctype列
kc.Phen$kc.eG4 <- factor(kc.Phen$kc.eG4,levels=c("eG4","non-eG4","other"))
my_comparisons = list(c("eG4","non-eG4"),c("non-eG4","other"),c("eG4","other"))

ggplot(data = kc.Phen,aes(x=kc.eG4,y=log2FoldChange,fill=kc.eG4)) +
  geom_violin(position = position_dodge(width = 1), scale = 'width') +
  geom_boxplot(position = position_dodge(width = 1), outlier.size = 0.6, width = 0.2, show.legend = FALSE) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(4.5,5,5.5),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(-3,6)) +
  scale_fill_manual(values = c("#A4A4A4","#3886B9","#D24E5D")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2](PhenDC3/DMSO))) +
  labs(title="Kc167 cell") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") 

tapply(kc.Phen$log2FoldChange,kc.Phen$kc.eG4, mean)
tapply(kc.Phen$log2FoldChange,kc.Phen$kc.eG4, median)

sample_data <- kc.Phen[kc.Phen$kc.eG4=="eG4",2]
# 单样本Wilcoxon检验(看大于0是否有显著性)
wilcox.test(sample_data, mu = 0,alternative = "greater") #p-value < 2.2e-16
sample_data <- kc.Phen[kc.Phen$kc.eG4=="non-eG4",2]
wilcox.test(sample_data, mu = 0,alternative = "greater") #p-value = 3.375e-05
sample_data <- kc.Phen[kc.Phen$kc.eG4=="other",2]
wilcox.test(sample_data, mu = 0,alternative = "less") #p-value < 2.2e-16 显著小于0

#### Kc启动子含有eG4\non-eG4\withouteG4差异倍数 ####
# sed '1d' /home/yuss/flyG4/data/ref/promoter2000.bed | bedtools intersect -a - -b /home/yuss/flyG4/result/PQS/001.2.dmel.pqs.bed -c > /home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter2000.pqs.bed
rm(list = ls());gc();rm(list = ls())
Num = "010.6."
kc.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.6.DE.kc.Phen.txt") %>% as.data.frame()
promoter.pqs <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter.pqs.bed") %>% as.data.frame()
kc.Phen$pqsnum <- promoter.pqs[match(kc.Phen$geneid,promoter.pqs$V4),6]
kc.Phen$pqs.type <- ifelse(kc.Phen$pqsnum==0,"other","pqs")
promoter.kcG4 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter2000.kc.bed") %>% as.data.frame()
kc.Phen$num <- promoter.kcG4[match(kc.Phen$geneid,promoter.kcG4$V4),6]
kc.Phen$kc.eG4 <- ifelse(kc.Phen$pqs.type=="pqs",ifelse(kc.Phen$num=="0","non-eG4","eG4"),"other") #增加kctype列
kc.Phen$kc.eG4 <- factor(kc.Phen$kc.eG4,levels=c("other","non-eG4","eG4"))
my_comparisons = list(c("eG4","non-eG4"),c("non-eG4","other"),c("eG4","other"))
ggplot(data = kc.Phen,aes(x=kc.eG4,y=log2FoldChange,fill=kc.eG4)) +
  geom_violin(position = position_dodge(width = 1), scale = 'width') +
  geom_boxplot(position = position_dodge(width = 1), outlier.size = 0.6, width = 0.2, show.legend = FALSE) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(4.5,5,5.5),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(-3,6)) +
  scale_fill_manual(values = c("#A4A4A4","#3886B9","#D24E5D")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2](PhenDC3/DMSO))) +
  labs(title="Kc167 cell (Promoter eG4)") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") 

tapply(kc.Phen$log2FoldChange,kc.Phen$kc.eG4, mean)
tapply(kc.Phen$log2FoldChange,kc.Phen$kc.eG4, median)

sample_data <- kc.Phen[kc.Phen$kc.eG4=="eG4",2]
# 单样本Wilcoxon检验(看大于0是否有显著性)
wilcox.test(sample_data, mu = 0,alternative = "greater") #p-value = 3.446e-14
sample_data <- kc.Phen[kc.Phen$kc.eG4=="non-eG4",2]
wilcox.test(sample_data, mu = 0,alternative = "greater") #p-value = 3.993e-16
sample_data <- kc.Phen[kc.Phen$kc.eG4=="other",2]
wilcox.test(sample_data, mu = 0,alternative = "less") #p-value = 6.713e-14 显著小于0

#### S2 ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "010.6."
s2.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.6.DE.s2.Phen.txt") %>% as.data.frame()
gene.pqs <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.2.gene.pqs.bed") %>% as.data.frame()
s2.Phen$pqsnum <- gene.pqs[match(s2.Phen$geneid,gene.pqs$V4),7]
s2.Phen$pqs.type <- ifelse(s2.Phen$pqsnum==0,"other","pqs")
gene.s2 <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.gene.s2.bed") %>% as.data.frame()
s2.Phen$num <- gene.s2[match(s2.Phen$geneid,gene.s2$V4),7]
s2.Phen$s2.eG4 <- ifelse(s2.Phen$pqs.type=="pqs",ifelse(s2.Phen$num=="0","non-eG4","eG4"),"other") #增加kctype列
s2.Phen$s2.eG4 <- factor(s2.Phen$s2.eG4,levels=c("eG4","non-eG4","other"))
my_comparisons = list(c("eG4","non-eG4"),c("non-eG4","other"),c("eG4","other"))

ggplot(data = s2.Phen,aes(x=s2.eG4,y=log2FoldChange,fill=s2.eG4)) +
  geom_violin(position = position_dodge(width = 1), scale = 'width') +
  geom_boxplot(position = position_dodge(width = 1), outlier.size = 0.6, width = 0.2, show.legend = FALSE) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(8.0,8.1,9.1),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(-5,10)) +
  scale_fill_manual(values = c("#A4A4A4","#3886B9","#D24E5D")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2](PhenDC3/DMSO))) +
  labs(title="S2 cell") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") 

tapply(s2.Phen$log2FoldChange,s2.Phen$s2.eG4, mean)
tapply(s2.Phen$log2FoldChange,s2.Phen$s2.eG4, median)

sample_data <- s2.Phen[s2.Phen$s2.eG4=="eG4",2]
# 单样本Wilcoxon检验(看大于0是否有显著性)
wilcox.test(sample_data, mu = 0,alternative = "greater") #p-value = 1.049e-13
sample_data <- s2.Phen[s2.Phen$s2.eG4=="non-eG4",2]
wilcox.test(sample_data, mu = 0,alternative = "greater") #p-value = 9.213e-14
sample_data <- s2.Phen[s2.Phen$s2.eG4=="other",2]
wilcox.test(sample_data, mu = 0,alternative = "greater") #p-value = 2.665e-05

#### S2启动子含有eG4\non-eG4\withouteG4差异倍数 ####
rm(list = ls());gc();rm(list = ls())
Num = "010.6."
s2.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.6.DE.s2.Phen.txt") %>% as.data.frame()
promoter.pqs <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter.pqs.bed") %>% as.data.frame()
s2.Phen$pqsnum <- promoter.pqs[match(s2.Phen$geneid,promoter.pqs$V4),6]
s2.Phen$pqs.type <- ifelse(s2.Phen$pqsnum==0,"other","pqs")
promoter.s2G4 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter2000.s2.bed") %>% as.data.frame()
s2.Phen$num <- promoter.s2G4[match(s2.Phen$geneid,promoter.s2G4$V4),6]
s2.Phen$s2.eG4 <- ifelse(s2.Phen$pqs.type=="pqs",ifelse(s2.Phen$num=="0","non-eG4","eG4"),"other") #增加type列
s2.Phen$s2.eG4 <- factor(s2.Phen$s2.eG4,levels=c("other","non-eG4","eG4"))
my_comparisons = list(c("eG4","non-eG4"),c("non-eG4","other"),c("eG4","other"))
ggplot(data = s2.Phen,aes(x=s2.eG4,y=log2FoldChange,fill=s2.eG4)) +
  geom_violin(position = position_dodge(width = 1), scale = 'width') +
  geom_boxplot(position = position_dodge(width = 1), outlier.size = 0.6, width = 0.2, show.legend = FALSE) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(6.5,7.9,8.4),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(-3,9)) +
  scale_fill_manual(values = c("#A4A4A4","#3886B9","#D24E5D")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2](PhenDC3/DMSO))) +
  labs(title="S2 cell (Promoter eG4)") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") 

# tapply(kc.Phen$log2FoldChange,kc.Phen$kc.eG4, mean)
# tapply(kc.Phen$log2FoldChange,kc.Phen$kc.eG4, median)
# 
# sample_data <- kc.Phen[kc.Phen$kc.eG4=="eG4",2]
# # 单样本Wilcoxon检验(看大于0是否有显著性)
# wilcox.test(sample_data, mu = 0,alternative = "greater") #p-value = 3.446e-14
# sample_data <- kc.Phen[kc.Phen$kc.eG4=="non-eG4",2]
# wilcox.test(sample_data, mu = 0,alternative = "greater") #p-value = 3.993e-16
# sample_data <- kc.Phen[kc.Phen$kc.eG4=="other",2]
# wilcox.test(sample_data, mu = 0,alternative = "less") #p-value = 6.713e-14 显著小于0
#*GO KEGG--------------------------------------------------------------------------
#https://mp.weixin.qq.com/s?__biz=Mzg5OTYzMzY5Ng==&mid=2247484957&idx=1&sn=f097a730785b1f78b9cf6648b95f041f&chksm=c0510152f7268844065326d18437e38795c4fc8615cc2bfb916e03a4c9708f334bb757df5019&scene=21#wechat_redirect
#### 1.Kc 上下调GOKEGG计算准备 ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "010.6."
library(clusterProfiler)
library(HDO.db)
library(DOSE)
library(dplyr)
library(data.table)

kc.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.6.DE.kc.Phen.txt") %>% as.data.frame()
df <- kc.Phen[kc.Phen$group=="Up",]
kcPhenup <- unlist(df$geneid)

##富集分析
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

enrich_kcPhenup <- enrich_go_kegg(kcPhenup)
enrich_kcPhenup_ego <- enrich_kcPhenup$ego_all
enrich_kcPhenup_kegg <- enrich_kcPhenup$kegg
enrich_kcPhenup_ego_ALL <- enrich_kcPhenup$ego_ALL
#### 2.Kc 上调GO可视化 ####
##可视化
dotplot(enrich_kcPhenup_ego_ALL,showCategory=10,split="ONTOLOGY",label_format=60,title="PhenDC3 up regulated(Kc) GO enrichment") + facet_grid(ONTOLOGY~., scale='free')
dotplot(enrich_kcPhenup_ego_ALL) + ggtitle("PhenDC3 up regulated(Kc) GO enrichment")
##showCategory指定展示的GO Terms的个数，默认展示显著富集的top10个
##label_format=60左边的名称一行显示60
##facet_grid(ONTOLOGY~., scale='free')将三个基因本体分开

##计算Enrichment Factor
eGo <- separate(data=enrich_kcPhenup_ego , col=GeneRatio,into = c("GR1", "GR2"), sep = "/") #劈分GeneRatio为2列（GR1、GR2）
eGo <- separate(data=eGo, col=BgRatio, into = c("BR1", "BR2"), sep = "/") #劈分BgRatio为2列（BR1、BR2）
eGo <- mutate(eGo, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) #计算Enrichment Factor 

##只用排名前10的term画图，先把BP、MF、CC三种类型各取top10，然后再组合在一起。
eGoBP <- eGo %>% 
  filter(ONTOLOGY=="BP") %>%
  filter(row_number() >= 1,row_number() <= 10)
eGoCC <- eGo %>% 
  filter(ONTOLOGY=="CC") %>%
  filter(row_number() >= 1,row_number() <= 10)
eGoMF <- eGo %>% 
  filter(ONTOLOGY=="MF") %>%
  filter(row_number() >= 1,row_number() <= 10)
##特定的一些GO_BP
eGoBP <- subset(eGo,eGo$Description%in%c("cation transport","cilium assembly","cilium organization",
"positive regulation of transcription by RNA polymerase II","cation transmembrane transport","metal ion transport",
"positive regulation of transcription, DNA-templated","circadian rhythm","adenylate cyclase-activating G protein-coupled receptor signaling pathway",
                                         "adrenergic receptor signaling pathway"
))
eGo10 <- rbind(eGoBP,eGoMF,eGoCC)
##画图
ggplot(eGo10,aes(enrichment_factor,Description)) + 
  geom_point(aes(size=Count,color=qvalue,shape=ONTOLOGY)) +
  scale_color_gradient(low="red",high = "green") + 
  labs(color="pvalue",size="Count", shape="Ontology",
       x="Enrichment Factor",y="GO term",title="GO enrichment") + 
  theme_bw() +facet_wrap( ~ ONTOLOGY) ##横向分屏

##可以用fct_reorder(factor(x), y, .fun = median, .desc = FALSE)函数（将x按照y的顺序排序）对绘图排序。
eGo10$WrappedDescription <- str_wrap(eGo10$Description, width = 60) 
ggplot(eGo10,aes(enrichment_factor, fct_reorder(factor(WrappedDescription), enrichment_factor))) + 
  geom_point(aes(size=Count,color=-1*log10(pvalue),shape=ONTOLOGY)) +
  scale_color_gradient(low = "purple", high = "yellow") + 
  labs(color=expression(-log[10](Pvalue)),size="Count", shape="Ontology",
       x="Enrichment Factor",y="",title="GO Enrichment of Upregulated\nGenes in Kc167 Cell") + 
  theme_bw() + facet_wrap( ~ ONTOLOGY,ncol= 1,scale='free') # + facet_wrap( ~ ONTOLOGY)
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"KcPhenUpGO.pdf"),
       device = "pdf",width = 7,height = 9)
# GO：0030707
# 卵巢卵泡细胞发育

# eGoBP$rounded_pvalue <- sprintf("%.2e", eGoBP$pvalue)
# eGoBP$Description <- factor(eGoBP$Description, levels = eGoBP$Description[order(eGoBP$pvalue)])
# ggplot(eGoBP,aes(Description, enrichment_factor,fill = pvalue))+
#   geom_bar(stat = "identity")+
#   geom_text(aes(label=rounded_pvalue, y=enrichment_factor+2),size=3)+
#   coord_flip() +
#   labs(x='',y='Enrichment Factor', title = 'GO enrichment of upregulated genes\n in Kc167 Cell (biological processes)')+
#   scale_fill_gradient(low = "#FF0000", high = "#00FF00") +
#   theme_bw()+
#   theme(panel.grid = element_blank(),
#         axis.ticks.y = element_blank(),
#         plot.title = element_text(hjust = 0.5, size = 10))
# 
# ggplot(eGoBP,aes(enrichment_factor, fct_reorder(factor(Description), enrichment_factor))) +
#   geom_point(aes(size=Count,color=-1*log10(pvalue),shape=ONTOLOGY)) +
#   scale_color_gradient(low="purple",high="yellow") +
#   labs(color=expression(-log[10](Pvalue)),size="Count", shape="Ontology",
#        x="Enrichment Factor",y="",title="GO enrichment of upregulated\ngenes in Kc167 Cell") +
#   theme_bw() + facet_wrap( ~ ONTOLOGY,ncol= 1,scale='free') 
# ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"KcPhenUpGO_BP.pdf"),
#        device = "pdf",width = 7,height = 4)

eGoBP$Description <- factor(eGoBP$Description, levels = eGoBP$Description[order(eGoBP$pvalue)])
eGoBP$WrappedDescription <- str_wrap(eGoBP$Description, width = 43)  # 设置每行最多显示15个字符
ggplot(eGoBP, aes(enrichment_factor, fct_reorder(factor(WrappedDescription), enrichment_factor))) +
  geom_point(aes(size = Count, color = -1 * log10(pvalue), shape = ONTOLOGY)) +
  scale_color_gradient(low = "purple", high = "yellow") +
  labs(color = expression(-log[10](Pvalue)), size = "Count", shape = "Ontology",
       x = "Enrichment Factor", y = "", title = "GO enrichment of upregulated\ngenes in Kc167 Cell") +
  theme_bw() 
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"KcPhenUpGO_BP.pdf"),
       device = "pdf",width = 5.8,height = 4)


#### 3.Kc 下调GO可视化 ####
##down
df <- kc.Phen[kc.Phen$group=="Down",]
kcPhendown <- unlist(df$geneid)
enrich_kcPhendown <- enrich_go_kegg(kcPhendown)
enrich_kcPhendown_ego <- enrich_kcPhendown$ego_all
enrich_kcPhendown_kegg <- enrich_kcPhendown$kegg
enrich_kcPhendown_ego_ALL <- enrich_kcPhendown$ego_ALL

##计算Enrichment Factor
eGo <- separate(data=enrich_kcPhendown_ego , col=GeneRatio,into = c("GR1", "GR2"), sep = "/") #劈分GeneRatio为2列（GR1、GR2）
eGo <- separate(data=eGo, col=BgRatio, into = c("BR1", "BR2"), sep = "/") #劈分BgRatio为2列（BR1、BR2）
eGo <- mutate(eGo, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) #计算Enrichment Factor 

##只用排名前10的term画图，先把BP、MF、CC三种类型各取top10，然后再组合在一起。
eGoBP <- eGo %>% 
  filter(ONTOLOGY=="BP") %>%
  filter(row_number() >= 1,row_number() <= 10)
eGoCC <- eGo %>% 
  filter(ONTOLOGY=="CC") %>%
  filter(row_number() >= 1,row_number() <= 10)
eGoMF <- eGo %>% 
  filter(ONTOLOGY=="MF") %>%
  filter(row_number() >= 1,row_number() <= 10)
##特定的一些GO_BP
eGoBP <- subset(eGo,eGo$Description%in%c("ubiquitin-dependent protein catabolic process","cellular protein catabolic process",
 "cellular response to unfolded protein","heat shock-mediated polytene chromosome puffing","polytene chromosome puffing",
 "protein folding","response to hypoxia","'de novo' post-translational protein folding","double-strand break repair via break-induced replication",
 "chaperone-mediated protein folding"
 ))
eGo10 <- rbind(eGoBP,eGoMF,eGoCC)
##画图
##可以用fct_reorder(factor(x), y, .fun = median, .desc = FALSE)函数（将x按照y的顺序排序）对绘图排序。
eGo10$WrappedDescription <- str_wrap(eGo10$Description, width = 60) 
ggplot(eGo10,aes(enrichment_factor, fct_reorder(factor(WrappedDescription), enrichment_factor))) + 
  geom_point(aes(size=Count,color=-1*log10(pvalue),shape=ONTOLOGY)) +
  scale_color_gradient(low = "purple", high = "yellow") + 
  labs(color=expression(-log[10](Pvalue)),size="Count", shape="Ontology",
       x="Enrichment Factor",y="",title="GO Enrichment of Downregulated\nGenes in Kc167 Cell") + 
  theme_bw() + facet_wrap( ~ ONTOLOGY,ncol= 1,scale='free') # + facet_wrap( ~ ONTOLOGY)
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"KcPhenDownGO.pdf"),
       device = "pdf",width = 7,height = 9)
##BP
eGoBP$Description <- factor(eGoBP$Description, levels = eGoBP$Description[order(eGoBP$pvalue)])
eGoBP$WrappedDescription <- str_wrap(eGoBP$Description, width = 43)  # 设置每行最多显示15个字符
ggplot(eGoBP, aes(enrichment_factor, fct_reorder(factor(WrappedDescription), enrichment_factor))) +
  geom_point(aes(size = Count, color = -1 * log10(pvalue), shape = ONTOLOGY)) +
  scale_color_gradient(low = "purple", high = "yellow") +
  labs(color = expression(-log[10](Pvalue)), size = "Count", shape = "Ontology",
       x = "Enrichment Factor", y = "", title = "GO enrichment of downregulated\ngenes in Kc167 Cell") +
  theme_bw() 
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"KcPhenDownGO_BP.pdf"),
       device = "pdf",width = 5.8,height = 4)
#### 上下调GO####
kceGoBP.down <- subset(eGo,eGo$Description%in%c("ubiquitin-dependent protein catabolic process","cellular protein catabolic process",
                                         "cellular response to unfolded protein","heat shock-mediated polytene chromosome puffing","polytene chromosome puffing",
                                         "protein folding","response to hypoxia","'de novo' post-translational protein folding","double-strand break repair via break-induced replication",
                                         "chaperone-mediated protein folding"
))%>%  
  mutate(Group = "PhenDC3 Down") %>%
  mutate(enrich = log10(pvalue))

kceGO <-rbind(kceGoBP.up,kceGoBP.down)
# eGoBP$Description <- factor(eGoBP$Description, levels = eGoBP$Description[order(eGoBP$pvalue)])
kceGO$Description=factor(kceGO$Description,levels = kceGO$Description[order(kceGO$enrichment_factor)])

ggplot(kceGO,aes(reorder(Description, enrichment_factor),enrichment_factor,fill=Group))+
  geom_col()+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color="black",size=10),
        axis.line.x = element_line(color='black'),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = 'none')+
  coord_flip()+
  geom_segment(aes(y=0, yend=0,x=0,xend=20.5))+
  geom_text(data = kceGO[which(kceGO$log10p>0),],aes(x=Description, y=-0.01, label=Description),
            hjust=1, size=4)+
  geom_text(data = kceGO[which(kceGO$log10p<0),],aes(x=Description, y=0.01, label=Description),
            hjust=0, size=4)+
  # geom_text(data = A[which(A$ratio>0),],aes(label=Padj),
  #           hjust=-0.1, size=4, color='red')+
  # geom_text(data = A[which(A$ratio<0),],aes(label=Padj),
  #           hjust=1.1, size=4, color="red")+
  scale_fill_manual(values = c("#1084A4",
                               "#8D4873"))+
  scale_x_discrete(expand = expansion(mult = c(0,0)))+
  ylim(-80, 50)+
  labs(x='', y='Ratio')

ggplot(kceGO,aes(x=log10p,y=Description,fill=Group))+
  geom_col(width = 0.7)+
  theme_classic()+
  theme(axis.title =element_text(size = 12,face = "plain", color = 'black'),axis.text =element_text(size = 12,face = "bold", color = 'black'))+	
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        legend.title = element_text(color="black",size=13,face = "plain"),
        legend.position = "top",
        legend.text = element_text(size=13,face = "plain"),
        axis.text = element_text(color="black",size=12,face = "plain"),
        axis.line.x = element_line(color='black'),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(size=14,face = "bold"))+
  ylab("")+
  geom_text(data = kceGO[which(kceGO$log10p>0),],aes(x=-0.01, y=Description, label=Description),
            hjust=1, size=4,color="black") +
  geom_text(data = kceGO[which(kceGO$log10p<0),],aes(x=0.01, y=Description, label=Description),
            hjust=0, size=4,color="black")+
  scale_fill_manual(values = c(
    'PhenDC3 Up'="#D6604D",
    'PhenDC3 Down'="#74ADD1"))+
  labs(x='-log(Pvalue)', y='')+
  ggtitle("Kc167 cell - PhenDC3 vs DMSO")
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"KcPhenUpDown.KEGG.pdf"),
       device = "pdf",width = 7.0,height = 5)

#### 4.Kc 上调KEGG可视化 ####
enrich_kcPhenup_kegg <- arrange(enrich_kcPhenup_kegg, enrich_kcPhenup_kegg$p.adjust)
enrich_kcPhenup_kegg$Description1 <- sapply(strsplit(enrich_kcPhenup_kegg$Description, " - "), function(x) x[1])
ggplot(enrich_kcPhenup_kegg[1:6,],aes(x=Count/87,y=Description1,colour=-1*log10(pvalue),size=Count))+
  geom_point()+
  scale_size(range = c(2, 8))+ #scale_size修改图中的点的大小，范围是2到8
  scale_color_gradient(low = "blue", high = "red")+
  theme_bw()+
  ylab("KEGG Pathway Terms")+
  xlab("Gene Ratio")+
  labs(color=expression(-log[10](PValue)), title = "PhenDC3 up regulated(Kc) kegg enrichment") +#expression函数改变样式，[]是用来添加下标，^是用来添加上标
  theme(text = element_text(size=15))


#### 5.Kc 下调KEGG可视化 ####
enrich_kcPhendown_kegg <- arrange(enrich_kcPhendown_kegg, enrich_kcPhendown_kegg$p.adjust)
enrich_kcPhendown_kegg$Description1 <- sapply(strsplit(enrich_kcPhendown_kegg$Description, " - "), function(x) x[1])
ggplot(enrich_kcPhendown_kegg[1:9,],aes(x=Count/71,y=Description1,colour=-1*log10(pvalue),size=Count))+
  geom_point()+
  scale_size(range = c(2, 8))+ #scale_size修改图中的点的大小，范围是2到8
  scale_color_gradient(low = "blue", high = "red")+
  theme_bw()+
  ylab("KEGG Pathway Terms")+
  xlab("Gene Ratio")+
  labs(color=expression(-log[10](PValue)), title = "PhenDC3 down regulated(Kc) kegg enrichment") +#expression函数改变样式，[]是用来添加下标，^是用来添加上标
  theme(text = element_text(size=15))

#### 6.Kc 上下调KEGG可视化 ####
kcupkegg <- enrich_kcPhenup_kegg[c(1,3:6), ] %>%
  mutate(Group = "PhenDC3 Up") %>%
  mutate(log10p = -log10(pvalue))

kcdownkegg <- enrich_kcPhendown_kegg[c(2:9), ] %>%
  mutate(Group = "PhenDC3 Down") %>%
  mutate(log10p = log10(pvalue))

kckegg <-rbind(kcupkegg,kcdownkegg)
kckegg$Description1=factor(kckegg$Description1,levels = kckegg$Description1[c(6:14,5:1)])

ggplot(kckegg,aes(x=log10p,y=Description1,fill=Group))+
  geom_col(width = 0.7)+
  theme_classic()+
  theme(axis.title =element_text(size = 12,face = "plain", color = 'black'),axis.text =element_text(size = 12,face = "bold", color = 'black'))+	
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        legend.title = element_text(color="black",size=13,face = "plain"),
        legend.position = "top",
        legend.text = element_text(size=13,face = "plain"),
        axis.text = element_text(color="black",size=12,face = "plain"),
        axis.line.x = element_line(color='black'),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(size=14,face = "bold"))+
  ylab("")+
  geom_text(data = kckegg[which(kckegg$log10p>0),],aes(x=-0.01, y=Description1, label=Description1),
            hjust=1, size=4,color="black") +
  geom_text(data = kckegg[which(kckegg$log10p<0),],aes(x=0.01, y=Description1, label=Description1),
            hjust=0, size=4,color="black")+
  scale_fill_manual(values = c(
    'PhenDC3 Up'="#D6604D",
    'PhenDC3 Down'="#74ADD1"))+
  labs(x='-log(Pvalue)', y='')+
  ggtitle("Kc167 cell - PhenDC3 vs DMSO")
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"KcPhenUpDown.KEGG.pdf"),
       device = "pdf",width = 7.0,height = 5)

#### 1.S2 上下调GOKEGG计算准备 ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "010.6."
library(clusterProfiler)
library(HDO.db)
library(DOSE)
library(dplyr)
library(data.table)

s2.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.6.DE.s2.Phen.txt") %>% as.data.frame()
df <- s2.Phen[s2.Phen$group=="Up",]
s2Phenup <- unlist(df$geneid)

##富集分析
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

enrich_s2Phenup1 <- enrich_go_kegg(s2Phenup)
enrich_s2Phenup_ego1 <- enrich_s2Phenup1$ego_all
enrich_s2Phenup_kegg1 <- enrich_s2Phenup1$kegg
enrich_s2Phenup_ego_ALL1 <- enrich_s2Phenup1$ego_ALL

#### 2.S2 上调GO可视化 ####
dotplot(enrich_s2Phenup_ego_ALL,showCategory=10,split="ONTOLOGY",label_format=60,title="PhenDC3 up regulated(S2) GO enrichment") + facet_grid(ONTOLOGY~., scale='free')
##计算Enrichment Factor
eGo <- separate(data=enrich_s2Phenup_ego , col=GeneRatio,into = c("GR1", "GR2"), sep = "/") #劈分GeneRatio为2列（GR1、GR2）
eGo <- separate(data=eGo, col=BgRatio, into = c("BR1", "BR2"), sep = "/") #劈分BgRatio为2列（BR1、BR2）
eGo <- mutate(eGo, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) #计算Enrichment Factor 

##只用排名前10的term画图，先把BP、MF、CC三种类型各取top10，然后再组合在一起。
eGoBP <- eGo %>% 
  filter(ONTOLOGY=="BP") %>%
  filter(row_number() >= 1,row_number() <= 10)
eGoCC <- eGo %>% 
  filter(ONTOLOGY=="CC") %>%
  filter(row_number() >= 1,row_number() <= 10)
eGoMF <- eGo %>% 
  filter(ONTOLOGY=="MF") %>%
  filter(row_number() >= 1,row_number() <= 10)
##特定的一些GO_BP
eGoBP <- subset(eGo,eGo$Description%in%c("cilium organization","cilium assembly","nucleosome assembly","nucleosome organization",
 "cation transport","inorganic ion transmembrane transport","chromatin remodeling","chromatin assembly","circadian rhythm",
 "DNA-templated transcription, initiation"
))
eGo10 <- rbind(eGoBP,eGoMF,eGoCC)

##可以用fct_reorder(factor(x), y, .fun = median, .desc = FALSE)函数（将x按照y的顺序排序）对绘图排序。
ggplot(eGo10,aes(enrichment_factor, fct_reorder(factor(Description), enrichment_factor))) + 
  geom_point(aes(size=Count,color=-1*log10(pvalue),shape=ONTOLOGY)) +
  scale_color_gradient(low = "purple", high = "yellow") + 
  labs(color=expression(-log[10](Pvalue)),size="Count", shape="Ontology",
       x="Enrichment Factor",y="",title="GO Enrichment of Upregulated\nGenes in S2 Cell") + 
  theme_bw() + facet_wrap( ~ ONTOLOGY,ncol= 1,scale='free') # + facet_wrap( ~ ONTOLOGY)
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"S2PhenUpGO.pdf"),
       device = "pdf",width = 7,height = 9)

##BP 根据enrichment_factor排列的
eGoBP$WrappedDescription <- str_wrap(eGoBP$Description, width = 43)  # 设置每行最多显示43个字符
ggplot(eGoBP, aes(enrichment_factor, fct_reorder(factor(WrappedDescription), enrichment_factor))) +
  geom_point(aes(size = Count, color = -1 * log10(pvalue), shape = ONTOLOGY)) +
  scale_color_gradient(low = "purple", high = "yellow") +
  labs(color = expression(-log[10](Pvalue)), size = "Count", shape = "Ontology",
       x = "Enrichment Factor", y = "", title = "GO enrichment of upregulated\ngenes in S2 cell") +
  theme_bw() 
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"S2PhenUpGO_BP.pdf"),
       device = "pdf",width = 5.8,height = 4)
#### 3.S2 下调GO可视化####
##down
df <- s2.Phen[s2.Phen$group=="Down",]
s2Phendown <- unlist(df$geneid)

##富集分析
enrich_s2Phendown <- enrich_go_kegg(s2Phendown)
enrich_s2Phendown_ego <- enrich_s2Phendown$ego_all
enrich_s2Phendown_kegg <- enrich_s2Phendown$kegg
enrich_s2Phendown_ego_ALL <- enrich_s2Phendown$ego_ALL
dotplot(enrich_s2Phendown_ego_ALL,showCategory=10,split="ONTOLOGY",label_format=60,title="PhenDC3 down regulated(S2) GO enrichment") + facet_grid(ONTOLOGY~., scale='free')

##计算Enrichment Factor
eGo <- separate(data=enrich_s2Phendown_ego , col=GeneRatio,into = c("GR1", "GR2"), sep = "/") #劈分GeneRatio为2列（GR1、GR2）
eGo <- separate(data=eGo, col=BgRatio, into = c("BR1", "BR2"), sep = "/") #劈分BgRatio为2列（BR1、BR2）
eGo <- mutate(eGo, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) #计算Enrichment Factor 

##只用排名前10的term画图，先把BP、MF、CC三种类型各取top10，然后再组合在一起。
eGoBP <- eGo %>% 
  filter(ONTOLOGY=="BP") %>%
  filter(row_number() >= 1,row_number() <= 10)
eGoCC <- eGo %>% 
  filter(ONTOLOGY=="CC") %>%
  filter(row_number() >= 1,row_number() <= 10)
eGoMF <- eGo %>% 
  filter(ONTOLOGY=="MF") %>%
  filter(row_number() >= 1,row_number() <= 10)
##特定的一些GO_BP
eGoBP <- subset(eGo,eGo$Description%in%c("cell adhesion","DNA replication","axon guidance","axonogenesis",
                                         "neuron recognition","regulation of fibroblast growth factor receptor signaling pathway",
                                         "immune system process","immune response","positive regulation of fibroblast growth factor receptor signaling pathway",
                                        "cell-cell adhesion"))
eGo10 <- rbind(eGoBP,eGoMF,eGoCC)

##可以用fct_reorder(factor(x), y, .fun = median, .desc = FALSE)函数（将x按照y的顺序排序）对绘图排序。
ggplot(eGo10,aes(enrichment_factor, fct_reorder(factor(Description), enrichment_factor))) + 
  geom_point(aes(size=Count,color=-1*log10(pvalue),shape=ONTOLOGY)) +
  scale_color_gradient(low = "purple", high = "yellow") + 
  labs(color=expression(-log[10](Pvalue)),size="Count", shape="Ontology",
       x="Enrichment Factor",y="",title="GO Enrichment of Downregulated\nGenes in S2 Cell") + 
  theme_bw() + facet_wrap( ~ ONTOLOGY,ncol= 1,scale='free') # + facet_wrap( ~ ONTOLOGY)
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"S2PhenDownGO.pdf"),
       device = "pdf",width = 7,height = 9)

eGoBP$WrappedDescription <- str_wrap(eGoBP$Description, width = 43)  # 设置每行最多显示43个字符
ggplot(eGoBP, aes(enrichment_factor, fct_reorder(factor(WrappedDescription), enrichment_factor))) +
  geom_point(aes(size = Count, color = -1 * log10(pvalue), shape = ONTOLOGY)) +
  scale_color_gradient(low = "purple", high = "yellow") +
  labs(color = expression(-log[10](Pvalue)), size = "Count", shape = "Ontology",
       x = "Enrichment Factor", y = "", title = "GO enrichment of downregulated\ngenes in S2 cell") +
  theme_bw() 
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"S2PhenDownGO_BP.pdf"),
       device = "pdf",width = 5.8,height = 4.4)
#### 4.S2 上调KEGG可视化####
enrich_s2Phenup_kegg <- arrange(enrich_s2Phenup_kegg, enrich_s2Phenup_kegg$p.adjust)
enrich_s2Phenup_kegg$Description1 <- sapply(strsplit(enrich_s2Phenup_kegg$Description, " - "), function(x) x[1])
ggplot(enrich_s2Phenup_kegg[1:10,],aes(x=Count/96,y=Description1,colour=-1*log10(pvalue),size=Count))+
  geom_point()+
  scale_size(range = c(2, 8))+ #scale_size修改图中的点的大小，范围是2到8
  scale_color_gradient(low = "blue", high = "red")+
  theme_bw()+
  ylab("KEGG Pathway Terms")+
  xlab("Gene Ratio")+
  labs(color=expression(-log[10](PValue)), title = "PhenDC3 up regulated(S2) kegg enrichment") +#expression函数改变样式，[]是用来添加下标，^是用来添加上标
  theme(text = element_text(size=15))

#### 5.S2 下调KEGG可视化####
enrich_s2Phendown_kegg <- arrange(enrich_s2Phendown_kegg, enrich_s2Phendown_kegg$p.adjust)
enrich_s2Phendown_kegg$Description1 <- sapply(strsplit(enrich_s2Phendown_kegg$Description, " - "), function(x) x[1])
ggplot(enrich_s2Phendown_kegg[1:12,],aes(x=Count/94,y=Description1,colour=-1*log10(pvalue),size=Count))+
  geom_point()+
  scale_size(range = c(2, 8))+ #scale_size修改图中的点的大小，范围是2到8
  scale_color_gradient(low = "blue", high = "red")+
  theme_bw()+
  ylab("KEGG Pathway Terms")+
  xlab("Gene Ratio")+
  labs(color=expression(-log[10](PValue)), title = "PhenDC3 down regulated(S2) kegg enrichment") +#expression函数改变样式，[]是用来添加下标，^是用来添加上标
  theme(text = element_text(size=15))

#### 6.S2 上下调KEGG可视化####  
s2upkegg <- enrich_s2Phenup_kegg[c(1,3:10), ] %>%
  mutate(Group = "PhenDC3 Up") %>%
  mutate(log10p = -log10(pvalue))

s2downkegg <-subset(enrich_s2Phendown_kegg,enrich_s2Phendown_kegg$Description1%in%c("Longevity regulating pathway",
  "Dorso-ventral axis formation","DNA replication","Wnt signaling pathway","ABC transporters","Hippo signaling pathway",                                                                                  
  "Endocytosis","Glycine, serine and threonine metabolism","Spliceosome","MAPK signaling pathway"))%>%
  mutate(Group = "PhenDC3 Down") %>%
  mutate(log10p = log10(pvalue))
s2downkegg <- s2downkegg[c(1:10),]
s2kegg <-rbind(s2upkegg,s2downkegg)
s2kegg$Description1=factor(s2kegg$Description1,levels = s2kegg$Description1[c(10:19,9:1)])

ggplot(s2kegg,aes(x=log10p,y=Description1,fill=Group))+
  geom_col(width = 0.7)+
  theme_classic()+
  theme(axis.title =element_text(size = 12,face = "plain", color = 'black'),axis.text =element_text(size = 12,face = "bold", color = 'black'))+	
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        legend.title = element_text(color="black",size=13,face = "plain"),
        legend.position = "top",
        legend.text = element_text(size=13,face = "plain"),
        axis.text = element_text(color="black",size=12,face = "plain"),
        axis.line.x = element_line(color='black'),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(size=14,face = "bold"))+
  ylab("")+
  geom_text(data = s2kegg[which(s2kegg$log10p>0),],aes(x=-0.01, y=Description1, label=Description1),
            hjust=1, size=4,color="black") +
  geom_text(data = s2kegg[which(s2kegg$log10p<0),],aes(x=0.01, y=Description1, label=Description1),
            hjust=0, size=4,color="black")+
  scale_fill_manual(values = c(
    'PhenDC3 Up'="#BF69A9",
    'PhenDC3 Down'="#8BC25F"))+
  labs(x='-log(Pvalue)', y='')+
  ggtitle("S2 cell - PhenDC3 vs DMSO")
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"S2PhenUpDown.KEGG.pdf"),
       device = "pdf",width = 8.2,height = 5)

#*--------------------------------------------------------------------
#*含有eG4的基因GO KEGG--------------------------------------------------------------------------
#### 1.Kc 上下调GOKEGG计算准备 ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "010.6."
library(clusterProfiler)
library(HDO.db)
library(DOSE)
library(dplyr)
library(data.table)

kc.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.6.DE.kc.Phen.txt") %>% as.data.frame()
gene.kc <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.gene.kc.bed") %>% as.data.frame()
kc.Phen$num <- gene.kc[match(kc.Phen$geneid,gene.kc$V4),7]
kc.Phen$kc.eG4 <- ifelse(kc.Phen$num==0,"no eG4","eG4")
df <- kc.Phen[kc.Phen$group=="Up",] ##655
df <- kc.Phen[kc.Phen$group=="Up"&kc.Phen$kc.eG4=="eG4",] ##300
kcPhenup <- unlist(df$geneid)

##富集分析
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

enrich_kcPhenup <- enrich_go_kegg(kcPhenup)
enrich_kcPhenup_ego <- enrich_kcPhenup$ego_all
enrich_kcPhenup_kegg <- enrich_kcPhenup$kegg
enrich_kcPhenup_ego_ALL <- enrich_kcPhenup$ego_ALL
#### 2.Kc 上调GO可视化 ####
##特定的一些GO_BP
eGoBP <- subset(enrich_kcPhenup_ego,enrich_kcPhenup_ego$Description%in%c("cell morphogenesis involved in neuron differentiation","positive regulation of transcription by RNA polymerase II",
                                         "positive regulation of transcription, DNA-templated","circadian rhythm","G protein-coupled receptor signaling pathway",
                                         "ovarian follicle cell development","autophagy","positive regulation of macromolecule biosynthetic process","synapse assembly","potassium ion transmembrane transport"
))
plot_Phen <- eGoBP

#### 3.Kc 下调GO可视化 ####
##down
df <- kc.Phen[kc.Phen$group=="Down",] #211
df <- kc.Phen[kc.Phen$group=="Down"&kc.Phen$kc.eG4=="eG4",] ##29
kcPhendown <- unlist(df$geneid)
enrich_kcPhendown <- enrich_go_kegg(kcPhendown)
enrich_kcPhendown_ego <- enrich_kcPhendown$ego_all
enrich_kcPhendown_kegg <- enrich_kcPhendown$kegg
enrich_kcPhendown_ego_ALL <- enrich_kcPhendown$ego_ALL
eGoBP <- subset(enrich_kcPhendown_ego,enrich_kcPhendown_ego$Description%in%c("regulation of DNA replication","activation of protein kinase activity","regulation of phosphorylation",
                                                                             "hemoglobin metabolic process","hemoglobin biosynthetic process","regulation of guanylate cyclase activity","cellular response to anoxia","regulation of mitotic cell cycle",
                                                                             "positive regulation of cellular protein metabolic process","ketone body catabolic process"
))
plot_DMSO <- eGoBP
#### 上下调GO ####
eGoBP_Phen_DMSO_plot<-rbind(plot_Phen,plot_DMSO)
eGoBP_Phen_DMSO_plot$WrappedDescription <- str_wrap(eGoBP_Phen_DMSO_plot$Description, width = 57)  # 设置每行最多显示15个字符
enrich_kcPhenup_ego$Group="PhenDC3"
enrich_kcPhendown_ego$Group="DMSO"
eGoBP_plot<-rbind(enrich_kcPhenup_ego,enrich_kcPhendown_ego)%>%filter(Description%in%c("cell morphogenesis involved in neuron differentiation","positive regulation of transcription by RNA polymerase II",
                                                                                       "positive regulation of transcription, DNA-templated","circadian rhythm","G protein-coupled receptor signaling pathway",
                                                                                       "ovarian follicle cell development","autophagy","positive regulation of macromolecule biosynthetic process","synapse assembly","potassium ion transmembrane transport",
                                                                                       "regulation of DNA replication","activation of protein kinase activity","regulation of phosphorylation",
                                                                                       "hemoglobin metabolic process","hemoglobin biosynthetic process","regulation of guanylate cyclase activity","cellular response to anoxia","regulation of mitotic cell cycle",
                                                                                       "positive regulation of cellular protein metabolic process","ketone body catabolic process"))
eGoBP_plot$Description=factor(eGoBP_plot$Description,levels = unique(eGoBP_Phen_DMSO_plot$Description)[20:1])
eGoBP_plot$WrappedDescription <- str_wrap(eGoBP_plot$Description, width = 57)  # 设置每行最多显示15个字符
eGoBP_plot$WrappedDescription=factor(eGoBP_plot$WrappedDescription,levels = unique(eGoBP_Phen_DMSO_plot$WrappedDescription)[20:1])
ggplot()+
  geom_point(data=eGoBP_plot[eGoBP_plot$pvalue<0.05,],aes(y=WrappedDescription,x=Group,color=pvalue,size=Count))+
  geom_point(data=eGoBP_plot[eGoBP_plot$pvalue>0.05,],aes(y=WrappedDescription,x=Group),color="grey",size=1)+
  scale_colour_gradient(low = "purple", high = "yellow", na.value = NA)+
  xlab("")+
  ylab("")+
  theme_minimal()+
  labs(color="Pvalue",size="Count")+
  theme(# 标题居中
    plot.title = element_text(hjust = 0.8,vjust = 1.5,size = 15,face = "bold"),
    axis.text = element_text(color = 'black',size=12,face = "plain"),
    axis.title.x  = element_text(color = 'black',size=12,face = "plain"),
    axis.title = element_text(color = 'black'),
    axis.title.y   = element_text(color = 'black',size=12,face = "plain"),
    legend.text = element_text(size=12,face = "plain"),
    legend.title = element_text(size=12,face = "plain"))+
  ggtitle("GO enrichment (BP) in Kc167 cell - PhenDC3 vs DMSO")
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"KcPhenUpDownGO_BP.pdf"),
       device = "pdf",width = 7.3,height = 6)

#### 4.Kc 上调KEGG可视化 ####
enrich_kcPhenup_kegg <- arrange(enrich_kcPhenup_kegg, enrich_kcPhenup_kegg$p.adjust)
enrich_kcPhenup_kegg$Description1 <- sapply(strsplit(enrich_kcPhenup_kegg$Description, " - "), function(x) x[1])

#### 5.Kc 下调KEGG可视化 ####
enrich_kcPhendown_kegg <- arrange(enrich_kcPhendown_kegg, enrich_kcPhendown_kegg$p.adjust)
enrich_kcPhendown_kegg$Description1 <- sapply(strsplit(enrich_kcPhendown_kegg$Description, " - "), function(x) x[1])

#### 6.Kc 上下调KEGG可视化 ####
kcupkegg <- enrich_kcPhenup_kegg[c(1,3:7), ] %>%
  mutate(Group = "PhenDC3 Up") %>%
  mutate(log10p = -log10(pvalue))

kcdownkegg <- enrich_kcPhendown_kegg[c(1:3), ] %>%
  mutate(Group = "PhenDC3 Down") %>%
  mutate(log10p = log10(pvalue))

kckegg <-rbind(kcupkegg,kcdownkegg)
kckegg$Description1[1] <- "Hippo signaling pathway - multiple species"
kckegg$Description1[5] <- "Hippo signaling pathway - fly"
kckegg$Description1=factor(kckegg$Description1,levels = kckegg$Description1[c(7:9,6:1)])

ggplot(kckegg,aes(x=log10p,y=Description1,fill=Group))+
  geom_col(width = 0.7)+
  theme_classic()+
  theme(axis.title =element_text(size = 12,face = "plain", color = 'black'),axis.text =element_text(size = 12,face = "bold", color = 'black'))+	
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        legend.title = element_text(color="black",size=13,face = "plain"),
        legend.position = "top",
        legend.text = element_text(size=13,face = "plain"),
        axis.text = element_text(color="black",size=12,face = "plain"),
        axis.line.x = element_line(color='black'),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(size=14,face = "bold"))+
  ylab("")+
  geom_text(data = kckegg[which(kckegg$log10p>0),],aes(x=-0.01, y=Description1, label=Description1),
            hjust=1, size=4,color="black") +
  geom_text(data = kckegg[which(kckegg$log10p<0),],aes(x=0.01, y=Description1, label=Description1),
            hjust=0, size=4,color="black")+
  scale_fill_manual(values = c(
    'PhenDC3 Up'="#D6604D",
    'PhenDC3 Down'="#74ADD1"))+
  labs(x='-log(Pvalue)', y='')+
  ggtitle("Kc167 cell - PhenDC3 vs DMSO")
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"eG4.KcPhenUpDown.KEGG.pdf"),
       device = "pdf",width = 7.5,height = 4.8)

#### 1.S2 上下调GOKEGG计算准备 ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "010.6."
library(clusterProfiler)
library(HDO.db)
library(DOSE)
library(dplyr)
library(data.table)

s2.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.6.DE.s2.Phen.txt") %>% as.data.frame()
gene.s2 <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.gene.s2.bed") %>% as.data.frame()
s2.Phen$num <- gene.s2[match(s2.Phen$geneid,gene.s2$V4),7]
s2.Phen$s2.eG4 <- ifelse(s2.Phen$num==0,"no eG4","eG4")
s2.Phen$s2.eG4 <- factor(s2.Phen$s2.eG4,levels = c("no eG4","eG4"))
df <- s2.Phen[s2.Phen$group=="Up",] ##1306
df <- s2.Phen[s2.Phen$group=="Up"&s2.Phen$s2.eG4=="eG4",] ##369
s2Phenup <- unlist(df$geneid)
##富集分析
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

enrich_s2Phenup <- enrich_go_kegg(s2Phenup)
enrich_s2Phenup_ego <- enrich_s2Phenup$ego_all
enrich_s2Phenup_kegg <- enrich_s2Phenup$kegg
enrich_s2Phenup_ego_ALL <- enrich_s2Phenup$ego_ALL
merged_df <- merge(enrich_kcPhenup_ego, enrich_s2Phenup_ego, by = "ID")

#### 2.S2 上调GO可视化 ####
##特定的一些GO_BP  
eGoBP <- subset(enrich_s2Phenup_ego,enrich_s2Phenup_ego$Description%in%c("negative regulation of transcription by RNA polymerase II","cell morphogenesis involved in differentiation","eye development",
                                                                         "male genitalia development","synapse assembly","male sex differentiation","positive regulation of transcription, DNA-templated","growth",
                                                                         "positive regulation of macromolecule biosynthetic process","potassium ion transmembrane transport"
))
plot_Phen <- eGoBP

#### 3.S2 下调GO可视化####
##down
df <- s2.Phen[s2.Phen$group=="Down",] ##740
df <- s2.Phen[s2.Phen$group=="Down"&s2.Phen$s2.eG4=="eG4",] ##253
s2Phendown <- unlist(df$geneid)

##富集分析
enrich_s2Phendown <- enrich_go_kegg(s2Phendown)
enrich_s2Phendown_ego <- enrich_s2Phendown$ego_all
enrich_s2Phendown_kegg <- enrich_s2Phendown$kegg
enrich_s2Phendown_ego_ALL <- enrich_s2Phendown$ego_ALL
eGoBP <- subset(enrich_s2Phendown_ego,enrich_s2Phendown_ego$Description%in%c("cell adhesion","axon guidance","sex differentiation","gonad development","positive regulation of cell cycle process",
                                                                             "negative regulation of cell differentiation",
                                         "regulation of fibroblast growth factor receptor signaling pathway","neuron recognition",
                                         "regulation of phosphorylation","immune response","regulation of DNA replication"
                                         ))
merged_df <- merge(enrich_kcPhendown_ego, enrich_s2Phendown_ego, by = "ID")
plot_DMSO <- eGoBP

eGoBP_Phen_DMSO_plot<-rbind(plot_Phen,plot_DMSO)
eGoBP_Phen_DMSO_plot$WrappedDescription <- str_wrap(eGoBP_Phen_DMSO_plot$Description, width = 57)  # 设置每行最多显示15个字符
enrich_s2Phenup_ego$Group="PhenDC3"
enrich_s2Phendown_ego$Group="DMSO"
eGoBP_plot<-rbind(enrich_s2Phenup_ego,enrich_s2Phendown_ego)%>%filter(Description%in%c("negative regulation of transcription by RNA polymerase II","cell morphogenesis involved in differentiation","eye development",
                                                                                       "male genitalia development","synapse assembly","male sex differentiation","positive regulation of transcription, DNA-templated","growth",
                                                                                       "positive regulation of macromolecule biosynthetic process","potassium ion transmembrane transport",
                                                                                       "cell adhesion","axon guidance","sex differentiation","gonad development","positive regulation of cell cycle process",
                                                                                       "negative regulation of cell differentiation",
                                                                                       "regulation of fibroblast growth factor receptor signaling pathway","neuron recognition",
                                                                                       "regulation of phosphorylation","immune response","regulation of DNA replication"))%>%filter(!(GeneRatio %in% c("41/217", "26/217", "24/217", "12/217")))
eGoBP_plot$WrappedDescription <- str_wrap(eGoBP_plot$Description, width = 57)  # 设置每行最多显示15个字符
eGoBP_plot$WrappedDescription=factor(eGoBP_plot$WrappedDescription,levels = unique(eGoBP_Phen_DMSO_plot$WrappedDescription)[20:1])
ggplot()+
  geom_point(data=eGoBP_plot[eGoBP_plot$pvalue<0.05,],aes(y=WrappedDescription,x=Group,color=pvalue,size=Count))+
  geom_point(data=eGoBP_plot[eGoBP_plot$pvalue>0.05,],aes(y=WrappedDescription,x=Group),color="grey",size=1)+
  scale_colour_gradient(low = "purple", high = "yellow", na.value = NA)+
  xlab("")+
  ylab("")+
  theme_minimal()+
  labs(color="Pvalue",size="Count")+
  theme(# 标题居中
    plot.title = element_text(hjust = 0.8,vjust = 1.5,size = 15,face = "bold"),
    axis.text = element_text(color = 'black',size=12,face = "plain"),
    axis.title.x  = element_text(color = 'black',size=12,face = "plain"),
    axis.title = element_text(color = 'black'),
    axis.title.y   = element_text(color = 'black',size=12,face = "plain"),
    legend.text = element_text(size=12,face = "plain"),
    legend.title = element_text(size=12,face = "plain"))+
  ggtitle("GO enrichment (BP) in S2 cell - PhenDC3 vs DMSO")
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"S2PhenUpDownGO_BP.pdf"),
       device = "pdf",width = 7.9,height = 6)


#### 4.S2 上调KEGG可视化####
enrich_s2Phenup_kegg <- arrange(enrich_s2Phenup_kegg, enrich_s2Phenup_kegg$p.adjust)
enrich_s2Phenup_kegg$Description1 <- sapply(strsplit(enrich_s2Phenup_kegg$Description, " - "), function(x) x[1])

#### 5.S2 下调KEGG可视化####
enrich_s2Phendown_kegg <- arrange(enrich_s2Phendown_kegg, enrich_s2Phendown_kegg$p.adjust)
enrich_s2Phendown_kegg$Description1 <- sapply(strsplit(enrich_s2Phendown_kegg$Description, " - "), function(x) x[1])


#### 6.S2 上下调KEGG可视化####  
s2upkegg <- enrich_s2Phenup_kegg[c(1:2,4:6), ] %>%
  mutate(Group = "PhenDC3 Up") %>%
  mutate(log10p = -log10(pvalue))
s2upkegg$Description1[1] <- "Hippo signaling pathway - multiple species"
s2upkegg$Description1[5] <- "Hippo signaling pathway - fly"
s2downkegg <- enrich_s2Phendown_kegg[c(1:9), ] %>%
  mutate(Group = "PhenDC3 Down") %>%
  mutate(log10p = log10(pvalue))

s2kegg <-rbind(s2upkegg,s2downkegg)

s2kegg$Description1=factor(s2kegg$Description1,s2kegg$Description1[c(6:14,5:1)])
ggplot(s2kegg,aes(x=log10p,y=Description1,fill=Group))+
  geom_col(width = 0.7)+
  theme_classic()+
  theme(axis.title =element_text(size = 12,face = "plain", color = 'black'),axis.text =element_text(size = 12,face = "bold", color = 'black'))+	
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        legend.title = element_text(color="black",size=13,face = "plain"),
        legend.position = "top",
        legend.text = element_text(size=13,face = "plain"),
        axis.text = element_text(color="black",size=12,face = "plain"),
        axis.line.x = element_line(color='black'),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(size=14,face = "bold"))+
  ylab("")+
  geom_text(data = s2kegg[which(s2kegg$log10p>0),],aes(x=-0.01, y=Description1, label=Description1),
            hjust=1, size=4,color="black") +
  geom_text(data = s2kegg[which(s2kegg$log10p<0),],aes(x=0.01, y=Description1, label=Description1),
            hjust=0, size=4,color="black")+
  scale_fill_manual(values = c(
    'PhenDC3 Up'="#BF69A9",
    'PhenDC3 Down'="#8BC25F"))+
  labs(x='-log(Pvalue)', y='')+
  ggtitle("S2 cell - PhenDC3 vs DMSO")
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"eG4.S2PhenUpDown.KEGG.pdf"),
       device = "pdf",width = 6.2,height = 5.5)

#*基因的启动子含有eG4的GO KEGG--------------------------------------------------------------------------
#### 1.Kc 上下调GOKEGG计算准备 ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "010.6."
library(clusterProfiler)
library(HDO.db)
library(DOSE)
library(dplyr)
library(data.table)

kc.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.6.DE.kc.Phen.txt") %>% as.data.frame()
promoter.kcG4 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter2000.kc.bed") %>% as.data.frame()
kc.Phen$num <- promoter.kcG4[match(kc.Phen$geneid,promoter.kcG4$V4),6]
kc.Phen$kc.eG4 <- ifelse(kc.Phen$num==0,"no eG4","eG4")
df <- kc.Phen[kc.Phen$group=="Up",] ##655
df <- kc.Phen[kc.Phen$group=="Up"&kc.Phen$kc.eG4=="eG4",] #227
kcPhenup <- unlist(df$geneid)

##富集分析
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

enrich_kcPhenup <- enrich_go_kegg(kcPhenup)
enrich_kcPhenup_ego <- enrich_kcPhenup$ego_all
enrich_kcPhenup_kegg <- enrich_kcPhenup$kegg
enrich_kcPhenup_ego_ALL <- enrich_kcPhenup$ego_ALL
#### 2.Kc 上调GO可视化 ####
##特定的一些GO_BP
eGoBP <- enrich_kcPhenup_ego[c(1:7,11:13),]
plot_Phen <- eGoBP
plot_Phen$Group="PhenDC3"
#### 3.Kc 下调GO可视化 ####
##down
df <- kc.Phen[kc.Phen$group=="Down",] #211
df <- kc.Phen[kc.Phen$group=="Down"&kc.Phen$kc.eG4=="eG4",] ##36
kcPhendown <- unlist(df$geneid)
enrich_kcPhendown <- enrich_go_kegg(kcPhendown)
enrich_kcPhendown_ego <- enrich_kcPhendown$ego_all
enrich_kcPhendown_kegg <- enrich_kcPhendown$kegg
enrich_kcPhendown_ego_ALL <- enrich_kcPhendown$ego_ALL
eGoBP <- enrich_kcPhendown_ego[1:9,]
plot_DMSO <- eGoBP
plot_DMSO$Group="DMSO"

eGoBP_Phen_DMSO_plot<-rbind(plot_Phen,plot_DMSO)

eGoBP_Phen_DMSO_plot$WrappedDescription <- str_wrap(eGoBP_Phen_DMSO_plot$Description, width = 57)  # 设置每行最多显示15个字符
eGoBP_Phen_DMSO_plot$WrappedDescription <- factor(eGoBP_Phen_DMSO_plot$WrappedDescription,levels = rev(eGoBP_Phen_DMSO_plot$WrappedDescription))

ggplot()+
  geom_point(data=eGoBP_Phen_DMSO_plot[eGoBP_Phen_DMSO_plot$pvalue<0.05,],aes(y=WrappedDescription,x=Group,color=pvalue,size=Count))+
  geom_point(data=eGoBP_Phen_DMSO_plot[eGoBP_Phen_DMSO_plot$pvalue>0.05,],aes(y=WrappedDescription,x=Group),color="grey",size=1)+
  scale_colour_gradient(low = "purple", high = "yellow", na.value = NA)+
  xlab("")+
  ylab("")+
  theme_minimal()+
  labs(color="Pvalue",size="Count")+
  theme(# 标题居中
    plot.title = element_text(hjust = 0.8,vjust = 1.5,size = 15,face = "bold"),
    axis.text = element_text(color = 'black',size=12,face = "plain"),
    axis.title.x  = element_text(color = 'black',size=12,face = "plain"),
    axis.title = element_text(color = 'black'),
    axis.title.y   = element_text(color = 'black',size=12,face = "plain"),
    legend.text = element_text(size=12,face = "plain"),
    legend.title = element_text(size=12,face = "plain"))+
  ggtitle("GO enrichment (BP) of PromotereG4DEG in Kc cell - PhenDC3 vs DMSO")
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"PromotereG4DEG.KcPhenUpDownGO_BP.pdf"),
       device = "pdf",width = 7.5,height = 6.5)

#### 1.S2 上下调GOKEGG计算准备 ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "010.6."
library(clusterProfiler)
library(HDO.db)
library(DOSE)
library(dplyr)
library(data.table)

s2.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.6.DE.s2.Phen.txt") %>% as.data.frame()
promoter.s2G4 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter2000.s2.bed") %>% as.data.frame()
s2.Phen$num <- promoter.s2G4[match(s2.Phen$geneid,promoter.s2G4$V4),6]
s2.Phen$s2.eG4 <- ifelse(s2.Phen$num==0,"no eG4","eG4")
df <- s2.Phen[s2.Phen$group=="Up",] ##1306
df <- s2.Phen[s2.Phen$group=="Up"&s2.Phen$s2.eG4=="eG4",] ##374
s2Phenup <- unlist(df$geneid)
##富集分析
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

enrich_s2Phenup <- enrich_go_kegg(s2Phenup)
enrich_s2Phenup_ego <- enrich_s2Phenup$ego_all ##没有BP
enrich_s2Phenup_kegg <- enrich_s2Phenup$kegg
enrich_s2Phenup_ego_ALL <- enrich_s2Phenup$ego_ALL
merged_df <- merge(enrich_kcPhenup_ego, enrich_s2Phenup_ego, by = "ID")

#### 2.S2 上调GO可视化 没有####


#### 3.S2 下调GO可视化####
##down
df <- s2.Phen[s2.Phen$group=="Down",] ##740
df <- s2.Phen[s2.Phen$group=="Down"&s2.Phen$s2.eG4=="eG4",] ##172
s2Phendown <- unlist(df$geneid)

##富集分析
enrich_s2Phendown <- enrich_go_kegg(s2Phendown)
enrich_s2Phendown_ego <- enrich_s2Phendown$ego_all
enrich_s2Phendown_kegg <- enrich_s2Phendown$kegg
enrich_s2Phendown_ego_ALL <- enrich_s2Phendown$ego_ALL
eGoBP <- subset(enrich_s2Phendown_ego,enrich_s2Phendown_ego$Description%in%c("cell adhesion","axon guidance","cell population proliferation","positive regulation of cell cycle",
                                                                             "regulation of fibroblast growth factor receptor signaling pathway","central nervous system development","neuron projection development","cell morphogenesis involved in differentiation",
                                                                             "compound eye morphogenesis","hemocyte proliferation"
))
merged_df <- merge(enrich_kcPhendown_ego, enrich_s2Phendown_ego, by = "ID") ##没有相同的
plot_DMSO <- eGoBP
plot_DMSO$WrappedDescription <- str_wrap(plot_DMSO$Description, width = 57)  # 设置每行最多显示15个字符
plot_DMSO$Group="DMSO"
#plot_DMSO$WrappedDescription <- factor(plot_DMSO$WrappedDescription, levels = plot_DMSO$WrappedDescription[order(plot_DMSO$pvalue)])
##逆序
plot_DMSO$WrappedDescription <- factor(plot_DMSO$WrappedDescription, levels = rev(plot_DMSO$WrappedDescription[order(plot_DMSO$pvalue)]))
ggplot()+
  geom_point(data=plot_DMSO[plot_DMSO$pvalue<0.05,],aes(y=WrappedDescription,x=Group,color=pvalue,size=Count))+
  scale_colour_gradient(low = "purple", high = "yellow", na.value = NA)+
  xlab("")+
  ylab("")+
  theme_minimal()+
  labs(color="Pvalue",size="Count")+
  theme(# 标题居中
    plot.title = element_text(hjust = 0.785,vjust = 1.5,size = 15,face = "bold"),
    axis.text = element_text(color = 'black',size=12,face = "plain"),
    axis.title.x  = element_text(color = 'black',size=12,face = "plain"),
    axis.title = element_text(color = 'black'),
    axis.title.y   = element_text(color = 'black',size=12,face = "plain"),
    legend.text = element_text(size=12,face = "plain"),
    legend.title = element_text(size=12,face = "plain"))+
  ggtitle("GO enrichment (BP) of PromotereG4DEG in S2 cell - PhenDC3 vs DMSO")
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"PromotereG4DEG.S2PhenUpDownGO_BP.pdf"),
       device = "pdf",width = 7.3,height = 4)

#*Kc更换基因背景集-----------------------------------------------------------------------
#### 1.kc筛选tpm>0的基因 ####
rm(list = ls());gc();rm(list = ls())
Num = "010.6."
kc.tpm <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.2.KcTpm.txt") %>% as.data.frame()
s2.tpm <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.2.S2Tpm.txt") %>% as.data.frame()
rownames(kc.tpm) <- kc.tpm$gene_id
kc.tpm <- kc.tpm[,c(12,1:3,7:8)]
kc.tpm$sum <- rowSums(kc.tpm[,2:6])
kc.tpm0gene <- kc.tpm[kc.tpm$sum>0,1]

#### 2.kc 含有eG4的上调基因 ####
#1.kc 含有eG4的上调基因
kc.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.6.DE.kc.Phen.txt") %>% as.data.frame()
gene.kc <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.gene.kc.bed") %>% as.data.frame()
kc.Phen$num <- gene.kc[match(kc.Phen$geneid,gene.kc$V4),7]
kc.Phen$kc.eG4 <- ifelse(kc.Phen$num==0,"no eG4","eG4")
df <- kc.Phen[kc.Phen$group=="Up",] ##655
df <- kc.Phen[kc.Phen$group=="Up"&kc.Phen$kc.eG4=="eG4",] ##300
kcPhenup <- unlist(df$geneid)

#2.转换ID 用ENTREZID做富集分析，比ENSEMBL id 做富集分析的显著性高一点点
kc.tpm0gene_trans <- bitr(kc.tpm0gene,fromType = 'ENSEMBL',toType = 'ENTREZID',OrgDb = "org.Dm.eg.db")  

#### 3.kc 含有eG4的上调基因 GO、KEGG ####
library(clusterProfiler)
library(HDO.db)
library(DOSE)
library(dplyr)
library(data.table)
library(org.Dm.eg.db)
enrich_go_kegg<-function(gene_id){
  result<-list()
  ensembl_2_entrezid<-bitr(gene_id,fromType ="ENSEMBL", toType = c("ENTREZID","SYMBOL"), OrgDb = org.Dm.eg.db)
  ego_ALL <- enrichGO(gene = as.character(ensembl_2_entrezid$ENTREZID),
                      OrgDb=org.Dm.eg.db,
                      keyType = "ENTREZID",
                      ont = "BP",
                      pAdjustMethod = "BH",
                      minGSSize = 1,
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05,
                      universe = kc.tpm0gene_trans$ENTREZID, #背景基因集
                      readable = TRUE)
  result$ego_all<-as.data.frame(ego_ALL)
  ensembl_2_kegg_id<-bitr_kegg(ensembl_2_entrezid$ENTREZID,fromType = "ncbi-geneid",toType = "kegg",organism = "dme")
  kegg<-enrichKEGG(gene = ensembl_2_entrezid$ENTREZID,organism = "dme", keyType = "ncbi-geneid",universe = kc.tpm0gene_trans$ENTREZID,pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
  result$kegg<-kegg@result
  result$ego_ALL<-ego_ALL
  return(result)
}
enrich_kcPhenup <- enrich_go_kegg(kcPhenup)
enrich_kcPhenup_ego <- enrich_kcPhenup$ego_all
enrich_kcPhenup_kegg <- enrich_kcPhenup$kegg
enrich_kcPhenup_ego_ALL <- enrich_kcPhenup$ego_ALL
#特定的通路
eGoBP <- subset(enrich_kcPhenup_ego,enrich_kcPhenup_ego$Description%in%c("cell morphogenesis involved in neuron differentiation","positive regulation of transcription by RNA polymerase II",
                                                                         "positive regulation of transcription, DNA-templated","circadian rhythm","G protein-coupled receptor signaling pathway",
                                                                         "ovarian follicle cell development","autophagy","positive regulation of macromolecule biosynthetic process","synapse assembly","potassium ion transmembrane transport"
))
plot_Phen <- eGoBP

#### 4.kc 含有eG4的下调基因 GO、KEGG ####
##1.含有eG4的下调基因
df <- kc.Phen[kc.Phen$group=="Down",] #211
df <- kc.Phen[kc.Phen$group=="Down"&kc.Phen$kc.eG4=="eG4",] ##29
kcPhendown <- unlist(df$geneid)

enrich_kcPhendown <- enrich_go_kegg(kcPhendown)
enrich_kcPhendown_ego <- enrich_kcPhendown$ego_all
enrich_kcPhendown_kegg <- enrich_kcPhendown$kegg
enrich_kcPhendown_ego_ALL <- enrich_kcPhendown$ego_ALL

eGoBP <- subset(enrich_kcPhendown_ego,enrich_kcPhendown_ego$Description%in%c("regulation of DNA replication","activation of protein kinase activity","regulation of phosphorylation",
                                                                             "hemoglobin metabolic process","hemoglobin biosynthetic process","regulation of guanylate cyclase activity","cellular response to anoxia","regulation of mitotic cell cycle",
                                                                             "positive regulation of cellular protein metabolic process","ketone body catabolic process"
))
plot_DMSO <- eGoBP

kcPhendown_trans <- bitr(kcPhendown,fromType = 'ENSEMBL',toType = 'ENTREZID',OrgDb = "org.Dm.eg.db")

# 分开做GO KEGG
# down_ego <- enrichGO(
#   gene  = kcPhendown_trans$ENTREZID,
#   keyType = "ENTREZID",
#   OrgDb   = org.Dm.eg.db,
#   universe = kc.tpm0gene_trans$ENTREZID, #背景基因集
#   ont     = "BP",
#   pAdjustMethod = "BH",
#   pvalueCutoff  = 0.05,
#   qvalueCutoff  = 0.05,
#   readable = TRUE,
#   minGSSize = 1)
# down_ego_df <- as.data.frame(down_ego)
# 
# kegg <- enrichKEGG(gene = kcPhendown_trans$ENTREZID,organism = 'dme',
#                    universe = kc.tpm0gene_trans$ENTREZID,
#                    keyType = "ncbi-geneid",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
# kegg_df <- kegg@result

#### 5.Kc 上下调GO、KEGG可视化 ####
#1.上下调GO 
eGoBP_Phen_DMSO_plot<-rbind(plot_Phen,plot_DMSO)
eGoBP_Phen_DMSO_plot$WrappedDescription <- str_wrap(eGoBP_Phen_DMSO_plot$Description, width = 57)  # 设置每行最多显示15个字符
enrich_kcPhenup_ego$Group="PhenDC3"
enrich_kcPhendown_ego$Group="DMSO"
eGoBP_plot<-rbind(enrich_kcPhenup_ego,enrich_kcPhendown_ego)%>%filter(Description%in%c("cell morphogenesis involved in neuron differentiation","positive regulation of transcription by RNA polymerase II",
                                                                                       "positive regulation of transcription, DNA-templated","circadian rhythm","G protein-coupled receptor signaling pathway",
                                                                                       "ovarian follicle cell development","autophagy","positive regulation of macromolecule biosynthetic process","synapse assembly","potassium ion transmembrane transport",
                                                                                       "regulation of DNA replication","activation of protein kinase activity","regulation of phosphorylation",
                                                                                       "hemoglobin metabolic process","hemoglobin biosynthetic process","regulation of guanylate cyclase activity","cellular response to anoxia","regulation of mitotic cell cycle",
                                                                                       "positive regulation of cellular protein metabolic process","ketone body catabolic process"))
eGoBP_plot$Description=factor(eGoBP_plot$Description,levels = unique(eGoBP_Phen_DMSO_plot$Description)[20:1])
eGoBP_plot$WrappedDescription <- str_wrap(eGoBP_plot$Description, width = 57)  # 设置每行最多显示15个字符
eGoBP_plot$WrappedDescription=factor(eGoBP_plot$WrappedDescription,levels = unique(eGoBP_Phen_DMSO_plot$WrappedDescription)[20:1])
ggplot()+
  geom_point(data=eGoBP_plot[eGoBP_plot$pvalue<0.05,],aes(y=WrappedDescription,x=Group,color=pvalue,size=Count))+
  geom_point(data=eGoBP_plot[eGoBP_plot$pvalue>0.05,],aes(y=WrappedDescription,x=Group),color="grey",size=1)+
  scale_colour_gradient(low = "purple", high = "yellow", na.value = NA)+
  xlab("")+
  ylab("")+
  theme_minimal()+
  labs(color="Pvalue",size="Count")+
  theme(# 标题居中
    plot.title = element_text(hjust = 0.8,vjust = 1.5,size = 15,face = "bold"),
    axis.text = element_text(color = 'black',size=12,face = "plain"),
    axis.title.x  = element_text(color = 'black',size=12,face = "plain"),
    axis.title = element_text(color = 'black'),
    axis.title.y   = element_text(color = 'black',size=12,face = "plain"),
    legend.text = element_text(size=12,face = "plain"),
    legend.title = element_text(size=12,face = "plain"))+
  ggtitle("GO enrichment (BP) in Kc167 cell - PhenDC3 vs DMSO")
ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"tpm0gene.KcPhenUpDownGO_BP.pdf"),
       device = "pdf",width = 7.3,height = 6)

#2.上下调KEGG
enrich_kcPhenup_kegg <- arrange(enrich_kcPhenup_kegg, enrich_kcPhenup_kegg$p.adjust)
enrich_kcPhenup_kegg$Description1 <- sapply(strsplit(enrich_kcPhenup_kegg$Description, " - "), function(x) x[1])
enrich_kcPhendown_kegg <- arrange(enrich_kcPhendown_kegg, enrich_kcPhendown_kegg$p.adjust)
enrich_kcPhendown_kegg$Description1 <- sapply(strsplit(enrich_kcPhendown_kegg$Description, " - "), function(x) x[1])

kcupkegg <- enrich_kcPhenup_kegg[c(1,3:7), ] %>%
  mutate(Group = "PhenDC3 Up") %>%
  mutate(log10p = -log10(pvalue))

kcdownkegg <- enrich_kcPhendown_kegg[c(1:3), ] %>%
  mutate(Group = "PhenDC3 Down") %>%
  mutate(log10p = log10(pvalue))

kckegg <-rbind(kcupkegg,kcdownkegg)
kckegg$Description1[1] <- "Hippo signaling pathway - multiple species"
kckegg$Description1[5] <- "Hippo signaling pathway - fly"
kckegg$Description1=factor(kckegg$Description1,levels = kckegg$Description1[c(7:9,6:1)])

ggplot(kckegg,aes(x=log10p,y=Description1,fill=Group))+
  geom_col(width = 0.7)+
  theme_classic()+
  theme(axis.title =element_text(size = 12,face = "plain", color = 'black'),axis.text =element_text(size = 12,face = "bold", color = 'black'))+	
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        legend.title = element_text(color="black",size=13,face = "plain"),
        legend.position = "top",
        legend.text = element_text(size=13,face = "plain"),
        axis.text = element_text(color="black",size=12,face = "plain"),
        axis.line.x = element_line(color='black'),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(size=14,face = "bold"))+
  ylab("")+
  geom_text(data = kckegg[which(kckegg$log10p>0),],aes(x=-0.01, y=Description1, label=Description1),
            hjust=1, size=4,color="black") +
  geom_text(data = kckegg[which(kckegg$log10p<0),],aes(x=0.01, y=Description1, label=Description1),
            hjust=0, size=4,color="black")+
  scale_fill_manual(values = c(
    'PhenDC3 Up'="#D6604D",
    'PhenDC3 Down'="#74ADD1"))+
  labs(x='-log(Pvalue)', y='')+
  ggtitle("Kc167 cell - PhenDC3 vs DMSO")

#*S2更换基因背景集-----------------------------------------------------------------------
#### 1.s2 筛选tpm>0的基因 ####
rm(list = ls());gc();rm(list = ls())
Num = "010.6."
s2.tpm <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.2.S2Tpm.txt") %>% as.data.frame()
rownames(s2.tpm) <- s2.tpm$gene_id
s2.tpm <- s2.tpm[,c(12,1:3,7:8)]
s2.tpm$sum <- rowSums(s2.tpm[,2:6])
s2.tpm0gene <- s2.tpm[s2.tpm$sum>0,1]

#### 2.s2 含有eG4的上调基因 ####
#1.s2 含有eG4的上调基因
s2.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.6.DE.s2.Phen.txt") %>% as.data.frame()
gene.s2 <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.gene.s2.bed") %>% as.data.frame()
s2.Phen$num <- gene.s2[match(s2.Phen$geneid,gene.s2$V4),7]
s2.Phen$s2.eG4 <- ifelse(s2.Phen$num==0,"no eG4","eG4")
s2.Phen$s2.eG4 <- factor(s2.Phen$s2.eG4,levels = c("no eG4","eG4"))
df <- s2.Phen[s2.Phen$group=="Up",] ##1306
df <- s2.Phen[s2.Phen$group=="Up"&s2.Phen$s2.eG4=="eG4",] ##369
s2Phenup <- unlist(df$geneid)

#2.转换ID 用ENTREZID做富集分析，比ENSEMBL id 做富集分析的显著性高一点点
s2.tpm0gene_trans <- bitr(s2.tpm0gene,fromType = 'ENSEMBL',toType = 'ENTREZID',OrgDb = "org.Dm.eg.db")  

#### 3.s2 含有eG4的上调基因 GO、KEGG ####
library(clusterProfiler)
library(HDO.db)
library(DOSE)
library(dplyr)
library(data.table)
library(org.Dm.eg.db)
enrich_go_kegg<-function(gene_id){
  result<-list()
  ensembl_2_entrezid<-bitr(gene_id,fromType ="ENSEMBL", toType = c("ENTREZID","SYMBOL"), OrgDb = org.Dm.eg.db)
  ego_ALL <- enrichGO(gene = as.character(ensembl_2_entrezid$ENTREZID),
                      OrgDb=org.Dm.eg.db,
                      keyType = "ENTREZID",
                      ont = "BP",
                      pAdjustMethod = "BH",
                      minGSSize = 1,
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05,
                      universe = s2.tpm0gene_trans$ENTREZID, #背景基因集
                      readable = TRUE)
  result$ego_all<-as.data.frame(ego_ALL)
  ensembl_2_kegg_id<-bitr_kegg(ensembl_2_entrezid$ENTREZID,fromType = "ncbi-geneid",toType = "kegg",organism = "dme")
  kegg<-enrichKEGG(gene = ensembl_2_entrezid$ENTREZID,organism = "dme", keyType = "ncbi-geneid",universe = s2.tpm0gene_trans$ENTREZID,pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
  result$kegg<-kegg@result
  result$ego_ALL<-ego_ALL
  return(result)
}
enrich_s2Phenup <- enrich_go_kegg(s2Phenup)
enrich_s2Phenup_ego <- enrich_s2Phenup$ego_all
enrich_s2Phenup_kegg <- enrich_s2Phenup$kegg
enrich_s2Phenup_ego_ALL <- enrich_s2Phenup$ego_ALL
#特定的通路
eGoBP <- subset(enrich_s2Phenup_ego,enrich_s2Phenup_ego$Description%in%c("negative regulation of transcription by RNA polymerase II","cell morphogenesis involved in differentiation","eye development",
                                                                         "male genitalia development","synapse assembly","male sex differentiation","positive regulation of transcription, DNA-templated","growth",
                                                                         "positive regulation of macromolecule biosynthetic process","potassium ion transmembrane transport"
))
plot_Phen <- eGoBP

#### 4.s2 含有eG4的下调基因 GO、KEGG ####
##1.含有eG4的下调基因
df <- s2.Phen[s2.Phen$group=="Down",] ##740
df <- s2.Phen[s2.Phen$group=="Down"&s2.Phen$s2.eG4=="eG4",] ##253
s2Phendown <- unlist(df$geneid)

##富集分析
enrich_s2Phendown <- enrich_go_kegg(s2Phendown)
enrich_s2Phendown_ego <- enrich_s2Phendown$ego_all
enrich_s2Phendown_kegg <- enrich_s2Phendown$kegg
enrich_s2Phendown_ego_ALL <- enrich_s2Phendown$ego_ALL
eGoBP <- subset(enrich_s2Phendown_ego,enrich_s2Phendown_ego$Description%in%c("cell adhesion","axon guidance","sex differentiation","gonad development","positive regulation of cell cycle process",
                                                                             "negative regulation of cell differentiation",
                                                                             "regulation of fibroblast growth factor receptor signaling pathway","neuron recognition",
                                                                             "regulation of phosphorylation","immune response","regulation of DNA replication"
))
plot_DMSO <- eGoBP

#### 5.s2 上下调GO、KEGG可视化 ####
#1.上下调GO 
eGoBP_Phen_DMSO_plot <- rbind(plot_Phen,plot_DMSO)
eGoBP_Phen_DMSO_plot$WrappedDescription <- str_wrap(eGoBP_Phen_DMSO_plot$Description, width = 57)  # 设置每行最多显示15个字符
enrich_s2Phenup_ego$Group="PhenDC3"
enrich_s2Phendown_ego$Group="DMSO"
eGoBP_plot<-rbind(enrich_s2Phenup_ego,enrich_s2Phendown_ego)%>%filter(Description%in%c("negative regulation of transcription by RNA polymerase II","cell morphogenesis involved in differentiation","eye development",
                                                                                       "male genitalia development","synapse assembly","male sex differentiation","positive regulation of transcription, DNA-templated","growth",
                                                                                       "positive regulation of macromolecule biosynthetic process","potassium ion transmembrane transport",
                                                                                       "cell adhesion","axon guidance","sex differentiation","gonad development","positive regulation of cell cycle process",
                                                                                       "negative regulation of cell differentiation",
                                                                                       "regulation of fibroblast growth factor receptor signaling pathway","neuron recognition",
                                                                                       "regulation of phosphorylation","immune response","regulation of DNA replication"))%>%filter(!(GeneRatio %in% c("41/217", "26/217", "24/217", "12/217")))
eGoBP_plot$WrappedDescription <- str_wrap(eGoBP_plot$Description, width = 57)  # 设置每行最多显示15个字符
eGoBP_plot$WrappedDescription=factor(eGoBP_plot$WrappedDescription,levels = unique(eGoBP_Phen_DMSO_plot$WrappedDescription)[20:1])
ggplot()+
  geom_point(data=eGoBP_plot[eGoBP_plot$pvalue<0.05,],aes(y=WrappedDescription,x=Group,color=pvalue,size=Count))+
  geom_point(data=eGoBP_plot[eGoBP_plot$pvalue>0.05,],aes(y=WrappedDescription,x=Group),color="grey",size=1)+
  scale_colour_gradient(low = "purple", high = "yellow", na.value = NA)+
  xlab("")+
  ylab("")+
  theme_minimal()+
  labs(color="Pvalue",size="Count")+
  theme(# 标题居中
    plot.title = element_text(hjust = 0.8,vjust = 1.5,size = 15,face = "bold"),
    axis.text = element_text(color = 'black',size=12,face = "plain"),
    axis.title.x  = element_text(color = 'black',size=12,face = "plain"),
    axis.title = element_text(color = 'black'),
    axis.title.y   = element_text(color = 'black',size=12,face = "plain"),
    legend.text = element_text(size=12,face = "plain"),
    legend.title = element_text(size=12,face = "plain"))+
  ggtitle("GO enrichment (BP) in S2 cell - PhenDC3 vs DMSO")

#2.上下调KEGG
enrich_s2Phenup_kegg <- arrange(enrich_s2Phenup_kegg, enrich_s2Phenup_kegg$p.adjust)
enrich_s2Phenup_kegg$Description1 <- sapply(strsplit(enrich_s2Phenup_kegg$Description, " - "), function(x) x[1])

enrich_s2Phendown_kegg <- arrange(enrich_s2Phendown_kegg, enrich_s2Phendown_kegg$p.adjust)
enrich_s2Phendown_kegg$Description1 <- sapply(strsplit(enrich_s2Phendown_kegg$Description, " - "), function(x) x[1])

s2upkegg <- enrich_s2Phenup_kegg[c(1:2,4), ] %>%
  mutate(Group = "PhenDC3 Up") %>%
  mutate(log10p = -log10(pvalue))
s2upkegg$Description1[1] <- "Hippo signaling pathway - multiple species"

s2downkegg <- enrich_s2Phendown_kegg[c(1:9), ] %>%
  mutate(Group = "PhenDC3 Down") %>%
  mutate(log10p = log10(pvalue))
s2downkegg$Description1[5] <- "Hippo signaling pathway - fly"
s2kegg <-rbind(s2upkegg,s2downkegg)

s2kegg$Description1=factor(s2kegg$Description1,s2kegg$Description1[c(4:12,3:1)])
ggplot(s2kegg,aes(x=log10p,y=Description1,fill=Group))+
  geom_col(width = 0.7)+
  theme_classic()+
  theme(axis.title =element_text(size = 12,face = "plain", color = 'black'),axis.text =element_text(size = 12,face = "bold", color = 'black'))+	
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        legend.title = element_text(color="black",size=13,face = "plain"),
        legend.position = "top",
        legend.text = element_text(size=13,face = "plain"),
        axis.text = element_text(color="black",size=12,face = "plain"),
        axis.line.x = element_line(color='black'),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(size=14,face = "bold"))+
  ylab("")+
  geom_text(data = s2kegg[which(s2kegg$log10p>0),],aes(x=-0.01, y=Description1, label=Description1),
            hjust=1, size=4,color="black") +
  geom_text(data = s2kegg[which(s2kegg$log10p<0),],aes(x=0.01, y=Description1, label=Description1),
            hjust=0, size=4,color="black")+
  scale_fill_manual(values = c(
    'PhenDC3 Up'="#BF69A9",
    'PhenDC3 Down'="#8BC25F"))+
  labs(x='-log(Pvalue)', y='')+
  ggtitle("S2 cell - PhenDC3 vs DMSO")
