rm(list = ls());gc();rm(list = ls())#清空
Num = "007.2."
#### 基因上含有不同eG4的表达水平 ####
##1.gene.bed与kc、s2细胞系中所有的G4取交集
##2.含有G4的基因和不含有G4的基因的表达水平
##脚本002.1.Daniel.gene.expression.R
#读取交集的表
gene.kc <- fread('/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.gene.kc.bed') %>% as.data.frame()
gene.s2 <- fread('/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.gene.s2.bed') %>% as.data.frame()
#合并两个表（染色体信息在前，交集信息在后）
library(dplyr)
gene.kc.s2 <- bind_cols(gene.kc,gene.s2$V7)
#修改列名
colnames(gene.kc.s2) <- c("chr","start","end","id","strand","gene_symbol","kc","s2")
#改$kc $s2
gene.kc.s2$kc <- ifelse(gene.kc.s2$kc==0,0,1)
gene.kc.s2$s2 <- ifelse(gene.kc.s2$s2==0,0,2)
#求kc和s2的sum,sum=1表示该gene有kc细胞系中的G4，sum=2表示该gene有s2细胞系中的G4
gene.kc.s2$sum <- rowSums(gene.kc.s2[,7:8])
write.table(gene.kc.s2,file = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/",Num,"gene.kc.s2.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)
gene.kcG4 <- gene.kc.s2[gene.kc.s2$sum==1,]
gene.s2G4 <- gene.kc.s2[gene.kc.s2$sum==2,]
gene.overlapG4 <- gene.kc.s2[gene.kc.s2$sum==3,]
gene.nonG4 <- gene.kc.s2[gene.kc.s2$sum==0,]
gene.mergeG4 <- gene.kc.s2[gene.kc.s2$sum!=0,]
#按行合并
gene.kcG4$type <- 'gene.kcG4'
gene.s2G4$type <- 'gene.s2G4'
gene.overlapG4$type <- 'gene.overlapG4'
gene.nonG4$type <- 'gene.nonG4'
gene.mergeG4$type <- 'gene.mergeG4'
gene.merge.all <- bind_rows(gene.overlapG4,gene.kcG4,gene.s2G4,gene.nonG4,gene.mergeG4)
#合并log2Foldchange
results <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.1.DE.kcvss2.txt") %>% as.data.frame()
gene.merge.all$logFC <- results[match(gene.merge.all$id,results$gene_id),3]
#RNAseq数据过滤这一步过滤了表达量低的reads,导致gene.merge.all与resault合并时，gene.merge.all中的logFC匹配不到，有NA值，这要去除掉
gene.merge.all=subset(gene.merge.all, logFC!= 'NA') #只去除logFC列含有NA值的行
gene.merge.all$type <- factor(gene.merge.all$type,levels = c("gene.nonG4","gene.kcG4","gene.s2G4","gene.overlapG4","gene.mergeG4"))
write.table(gene.merge.all,file = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/",Num,"gene.merge.all.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)
#画图
my_comparisons = list(c("gene.nonG4","gene.kcG4"),c("gene.kcG4","gene.s2G4"),
                      c("gene.s2G4","gene.overlapG4"),c("gene.overlapG4","gene.mergeG4"))
b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c('#4d4d4d',b)
median.logFC <- aggregate(logFC ~  type, gene.merge.all, median)
library(ggpubr)
ggplot(data = gene.merge.all,aes(x=type,y=logFC,color=type)) +
  geom_boxplot() +
  geom_point()+
  geom_jitter(width = 0.25)+ #点是抖动的且没有超出箱线图的宽度
  scale_color_manual(values = color)+ #更改箱线图箱子的颜色
  theme_bw() +
  stat_compare_means(comparisons = my_comparisons,label.y = c(rep(c(14.5,17),3)),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4) +
  ylab(bquote(Log[2]~Fold~Change))+
  labs(title="Kc167 vs S2") +
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") + #没有图例
  coord_cartesian(ylim = c(-16,19)) +
  scale_x_discrete(labels = c("non-eG4","kcspecific","s2specific","overlap","merge"))
ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"GeneExpression_With_eG4.pdf"),
       device = "pdf",width = 7,height = 5.8)  

gene.kcallG4 <- gene.kc.s2[gene.kc.s2$kc==1,]
gene.kcallG4$type <- "kcall"
gene.kcallnonG4 <- gene.kc.s2[gene.kc.s2$kc==0,]
gene.kcallnonG4$type <- "nonkcG4"
gene.s2allG4 <- gene.kc.s2[gene.kc.s2$s2==2,]
gene.s2allG4$type <- "s2all"
gene.s2allnonG4 <- gene.kc.s2[gene.kc.s2$s2==0,]
gene.s2allnonG4$type <- "nons2G4"
df <- bind_rows(gene.kcallG4,gene.s2allG4,gene.kcallnonG4,gene.s2allnonG4)
df$logFC <- results[match(df$id,results$gene_id),3]
df <- subset(df, logFC!='NA')
df$type <- factor(df$type,levels = c("nonkcG4","kcall","nons2G4","s2all"))
my_comparisons = list(c("nonkcG4","kcall"),c("nons2G4","s2all"))
b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c('#4d4d4d',b)
# median.logFC <- aggregate(logFC ~  type, gene.merge.all, median)
library(ggpubr)
ggplot(data = df,aes(x=type,y=logFC,color=type)) +
  geom_boxplot() +
  geom_point()+
  geom_jitter(width = 0.25)+ #点是抖动的且没有超出箱线图的宽度
  scale_color_manual(values = color)+ #更改箱线图箱子的颜色
  theme_bw() +
  stat_compare_means(comparisons = my_comparisons,label.y = c(rep(14.5,17),2),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4) +
  ylab(bquote(Log[2]~Fold~Change))+
  labs(title="Kc167 vs S2") +
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") + #没有图例
  coord_cartesian(ylim = c(-16,19)) +
  scale_x_discrete(labels = c("non-eG4","kcalleG4","non-eG4","s2alleG4"))


#### X染色体区域基因上含有不同eG4的表达水平 ####
my_comparisons = list(c("gene.nonG4","gene.kcG4"),c("gene.kcG4","gene.s2G4"),
                      c("gene.s2G4","gene.overlapG4"),c("gene.overlapG4","gene.mergeG4"))
gene.merge.all.X <- gene.merge.all[gene.merge.all$chr=="X",]
ggplot(data = gene.merge.all.X,aes(x=type,y=logFC,color=type)) +
  geom_boxplot() +
  geom_point()+
  geom_jitter(width = 0.25)+ #点是抖动的且没有超出箱线图的宽度
  scale_color_manual(values = color)+ #更改箱线图箱子的颜色
  theme_bw() +
  stat_compare_means(comparisons = my_comparisons,label.y = c(rep(c(14.5,17),2)),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4) +
  ylab(bquote(Log[2]~Fold~Change))+
  labs(title="Kc167 vs S2(X chromosome)") +
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") + #没有图例
  coord_cartesian(ylim = c(-16,19)) +
  scale_x_discrete(labels = c("non-eG4","kcspecific","s2specific","overlap","merge"))

df.X <- df[df$chr=="X",]
my_comparisons = list(c("nonkcG4","kcall"),c("nons2G4","s2all"))
# my_comparisons = list(c("gene.nonG4","kcall"),c("kcall","s2all"))
ggplot(data = df.X,aes(x=type,y=logFC,color=type)) +
  geom_boxplot() +
  geom_point()+
  geom_jitter(width = 0.25)+ #点是抖动的且没有超出箱线图的宽度
  scale_color_manual(values = color)+ #更改箱线图箱子的颜色
  theme_bw() +
  stat_compare_means(comparisons = my_comparisons,label.y = c(14.5,17),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4) +
  ylab(bquote(Log[2]~Fold~Change))+
  labs(title="Kc167 vs S2(X chromosome)") +
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") + #没有图例
  coord_cartesian(ylim = c(-16,19)) +
 scale_x_discrete(labels = c("non-eG4","kcalleG4","non-eG4","s2alleG4"))

#### 在promoter基因上含有不同eG4的表达水平 #### 
rm(list = ls());gc();rm(list = ls())#清空
Num = "007.2."

# Linux
# gene.sort <- fread('/home/yuss/flyG4/data/ref/gene.sort.bed')
# gene.sort.p <- gene.sort[gene.sort$V5=='+',] ##提取正链的信息
# gene.sort.m <- gene.sort[gene.sort$V5=='-',]
# promoter.p <- data.frame(gene.sort.p$V1,gene.sort.p$V2-1500,gene.sort.p$V2+1500,gene.sort.p$V4,gene.sort.p$V5)
# colnames(promoter.p) <- c("chr","start","end","gene_id","strand")
# promoter.m <- data.frame(gene.sort.m$V1,gene.sort.m$V3-1500,gene.sort.m$V3+1500,gene.sort.m$V4,gene.sort.m$V5)
# colnames(promoter.m) <- c("chr","start","end","gene_id","strand")
# promoter <- rbind(promoter.p,promoter.m)
# promoter$start <- ifelse(promoter$start <0, 0, promoter$start)
# write.table(promoter,file = '/home/yuss/flyG4/data/ref/promoter.bed',
#             sep = '\t',col.names = T,row.names = F,quote = F)
# ##2.promoter.bed与kc、s2细胞系中所有的G4取交集，若有交集表示启动子上的该基因片段有G4
# sed '1d' /home/yuss/flyG4/result/PQS/001.2.kc_all.bed > /home/yuss/flyG4/data/ref/kc_all.bed ##删除信息第一行
# sed '1d' /home/yuss/flyG4/result/PQS/001.2.s2_all.bed > /home/yuss/flyG4/data/ref/s2_all.bed
# sed '1d' /home/yuss/flyG4/data/ref/promoter.bed | bedtools intersect -a - -b /home/yuss/flyG4/data/ref/kc_all.bed -c > /home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.promoter.kc.bed
# sed '1d' /home/yuss/flyG4/data/ref/promoter.bed | bedtools intersect -a - -b /home/yuss/flyG4/data/ref/s2_all.bed -c > /home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.promoter.s2.bed

##3.读取promoter与s2,kc G4交集的文件
promoter.kc <- fread('/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.promoter.kc.bed') %>% as.data.frame()
promoter.s2 <- fread('/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.promoter.s2.bed') %>% as.data.frame()
#合并两个表
promoter.kc.s2 <- bind_cols(promoter.kc,promoter.s2$V6)
colnames(promoter.kc.s2) <- c("chr","start","end","gene_id","strand","kc","s2")
promoter.kc.s2$kc <- ifelse(promoter.kc.s2$kc==0,0,1)
promoter.kc.s2$s2 <- ifelse(promoter.kc.s2$s2==0,0,2)
#求sum
promoter.kc.s2$sum <- rowSums(promoter.kc.s2[,6:7])
write.table(promoter.kc.s2,file = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/",Num,"promoter.kc.s2.txt"),
sep = '\t',col.names = T,row.names = F,quote = F)
promoter.kcG4 <- promoter.kc.s2[promoter.kc.s2$sum==1,]
promoter.s2G4 <- promoter.kc.s2[promoter.kc.s2$sum==2,]
promoter.overlapG4 <- promoter.kc.s2[promoter.kc.s2$sum==3,]
promoter.non_G4 <- promoter.kc.s2[promoter.kc.s2$sum==0,]
promoter.mergeG4 <- promoter.kc.s2[promoter.kc.s2$sum!=0,]
#按行合并
promoter.kcG4$type <- "promoter.kc_eG4"
promoter.s2G4$type <- "promoter.s2_eG4"
promoter.overlapG4$type <- "promoter.overlap_eG4"
promoter.mergeG4$type <- "promoter.merge_eG4"
promoter.non_G4$type <- "promoter.non_eG4"
promoter.merge.all <- bind_rows(promoter.kcG4,promoter.s2G4,promoter.overlapG4,promoter.mergeG4,promoter.non_G4,promoter.non_G4)
#合并表达量
results <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.1.DE.kcvss2.txt") %>% as.data.frame()
promoter.merge.all$log2FoldChange <- results[match(promoter.merge.all$gene_id,results$gene_id),3]
promoter.merge.all <- subset(promoter.merge.all,log2FoldChange!='NA')
promoter.merge.all$type <- factor(promoter.merge.all$type,levels = c("promoter.non_eG4","promoter.kc_eG4",
                                                                     "promoter.s2_eG4","promoter.overlap_eG4",
                                                                     "promoter.merge_eG4"))
write.table(promoter.merge.all,file = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/",Num,"promoter.merge.all.txt"),
sep = '\t',col.names = T,row.names = F,quote = F)

my_comparisons = list(c("promoter.non_eG4","promoter.kc_eG4"),c("promoter.kc_eG4","promoter.s2_eG4"),
                      c("promoter.s2_eG4","promoter.overlap_eG4"),c("promoter.overlap_eG4","promoter.merge_eG4"))
b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(5)
color <- c('#4d4d4d',b)

library(ggpubr)
ggplot(data = promoter.merge.all,aes(x=type,y=log2FoldChange,color=type)) +
  geom_boxplot() +
  geom_point()+
  geom_jitter(width = 0.25)+ #点是抖动的且没有超出箱线图的宽度
  scale_color_manual(values = color)+ #更改箱线图箱子的颜色
  theme_bw() +
  stat_compare_means(comparisons = my_comparisons,label.y = c(rep(c(14.5,17),2)),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4) +
  ylab(bquote(Log[2]~Fold~Change))+
  labs(title="Promoter genes (Kc167 vs S2)") +
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") + #没有图例
  coord_cartesian(ylim = c(-16,19)) +
  scale_x_discrete(labels = c("non-eG4","kcspecific","s2specific","overlap","merge"))
#ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"PromoterGeneExpression_With_eG4.pdf"),
       device = "pdf",width = 7,height = 5.8) 

#### X染色体区域在promoter基因上含有不同eG4的表达水平 #### 
promoter.merge.all.X <- promoter.merge.all[promoter.merge.all$chr=="X",]
ggplot(data = promoter.merge.all.X,aes(x=type,y=log2FoldChange,color=type)) +
  geom_boxplot() +
  geom_point()+
  geom_jitter(width = 0.25)+ #点是抖动的且没有超出箱线图的宽度
  scale_color_manual(values = color)+ #更改箱线图箱子的颜色
  theme_bw() +
  stat_compare_means(comparisons = my_comparisons,label.y = c(rep(c(14.5,17),2)),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4) +
  ylab(bquote(Log[2]~Fold~Change))+
  labs(title="Promoter genes at X chromosome(Kc167 vs S2)") +
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") + #没有图例

  coord_cartesian(ylim = c(-16,19)) +
  scale_x_discrete(labels = c("non-eG4","kcspecific","s2specific","overlap","merge"))

promoter.kcall <- promoter.kc.s2[promoter.kc.s2$kc==1,]
promoter.kcall$type <- "kcalleG4"
promoter.kcnonG4 <- promoter.kc.s2[promoter.kc.s2$kc==0,]
promoter.kcnonG4$type <- "kcnoneG4"
promoter.s2all <- promoter.kc.s2[promoter.kc.s2$s2==2,]
promoter.s2all$type <- "s2alleG4"
promoter.s2nonG4 <- promoter.kc.s2[promoter.kc.s2$s2==0,]
promoter.s2nonG4$type <- "s2noneG4"

df <- bind_rows(promoter.kcall,promoter.kcnonG4,promoter.s2all,promoter.s2nonG4)
results <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.1.DE.kcvss2.txt") %>% as.data.frame()
df$logFC <- results[match(df$gene_id,results$gene_id),3]
df <- subset(df,logFC!="NA")
df$type <- factor(df$type,levels = c("kcnoneG4","kcalleG4","s2noneG4","s2alleG4"))
my_comparisons <- list(c("kcnoneG4","kcalleG4"),c("s2noneG4","s2alleG4"))
ggplot(data = df,aes(x=type,y=logFC,color=type)) +
  geom_boxplot() +
  geom_point()+
  geom_jitter(width = 0.25)+ #点是抖动的且没有超出箱线图的宽度
  scale_color_manual(values = color)+ #更改箱线图箱子的颜色
  theme_bw() +
  stat_compare_means(comparisons = my_comparisons,label.y = c(rep(c(14.5,17),2)),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4) +
  ylab(bquote(Log[2]~Fold~Change))+
  labs(title="Promoter genes (Kc167 vs S2)") +
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") + #没有图例
  coord_cartesian(ylim = c(-16,19))# +
  scale_x_discrete(labels = c("non-eG4","kcspecific","s2specific","overlap","merge"))

  
df.X <- df[df$chr=="X",]
ggplot(data = df.X,aes(x=type,y=logFC,color=type)) +
  geom_boxplot() +
  geom_point()+
  geom_jitter(width = 0.25)+ #点是抖动的且没有超出箱线图的宽度
  scale_color_manual(values = color)+ #更改箱线图箱子的颜色
  theme_bw() +
  stat_compare_means(comparisons = my_comparisons,label.y = c(rep(c(14.5,17),2)),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4) +
  ylab(bquote(Log[2]~Fold~Change))+
  labs(title="Promoter genes at X chromosomes(Kc167 vs S2)") +
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") + #没有图例
  coord_cartesian(ylim = c(-16,19))# +
scale_x_discrete(labels = c("non-eG4","kcspecific","s2specific","overlap","merge"))