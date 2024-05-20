#*在上下调表达的基因中含有G4的基因的数量---------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())#清空
Num = "007.2."
gene.merge.all <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.2.gene.merge.all.txt") %>% as.data.frame()
results <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.1.DE.kcvss2.txt") %>% as.data.frame()
gene.merge.all$signifi <- results[match(gene.merge.all$id,results$gene_id),10]
table(gene.merge.all$type,gene.merge.all$signifi)
gene.merge.all.signifi <- table(gene.merge.all$type,gene.merge.all$signifi) %>% as.data.frame()
library(tidyr)
gene.merge.all.signifi1 <- spread(gene.merge.all.signifi,Var2,Freq)
colnames(gene.merge.all.signifi1)[1] <- "type"
gene.merge.all.signifi1$type <- factor(gene.merge.all.signifi1$type,levels = c("gene.nonG4","gene.kcG4","gene.s2G4","gene.kcs2G4"))
b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c('#4d4d4d',b)
##在Kc细胞中上调表达的基因含有eG4的数量
library(ggplot2)
ggplot(gene.merge.all.signifi1,aes(x=type,y=Up,fill=type)) +
  geom_bar(stat = 'identity',position = position_dodge(0.7),width = 0.7,color ='white') +
  coord_cartesian(ylim = c(0,1800)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = color) +
  geom_text(aes(label=Up,vjust = -0.5),color="black",size=4) +
  cowplot::theme_half_open() +
  ylab("Number of eG4s genes in Kc167 cells") + ##通过bquote函数给图标签添加上下标
  theme(axis.title.y = element_text(size = 16), ##y坐标轴标题字体大小
        axis.title.x = element_blank(), ##删除x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.position = "none") +
  scale_x_discrete(labels = c("gene.non_eG4","gene.kc_eG4","gene.s2_eG4","gene.overlap_eG4"))

##在S2细胞中上调表达的基因含有eG4的数量
library(ggplot2)
ggplot(gene.merge.all.signifi1,aes(x=type,y=Down,fill=type)) +
  geom_bar(stat = 'identity',position = position_dodge(0.7),width = 0.7,color ='white') +
  coord_cartesian(ylim = c(0,1900)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = color) +
  geom_text(aes(label=Down,vjust = -0.5),color="black",size=4) +
  cowplot::theme_half_open() +
  ylab("Number of eG4s genes in S2 cells") + ##通过bquote函数给图标签添加上下标
  theme(axis.title.y = element_text(size = 16), ##y坐标轴标题字体大小
        axis.title.x = element_blank(), ##删除x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.position = "none") +
  scale_x_discrete(labels = c("gene.non_eG4","gene.kc_eG4","gene.s2_eG4","gene.overlap_eG4"))

#*在上下调表达的基因中在X染色体上含有G4的基因的数量-----------------------------------------------------------
rm(list = ls());gc();rm(list = ls())#清空
Num = "007.2."
gene.merge.all <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.2.gene.merge.all.txt") %>% as.data.frame()
results <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.1.DE.kcvss2.txt") %>% as.data.frame()
gene.merge.all.X <- gene.merge.all[gene.merge.all$chr=="X",]
gene.merge.all.X$signifi <- results[match(gene.merge.all.X$id,results$gene_id),10]
table(gene.merge.all.X$type,gene.merge.all.X$signifi)
