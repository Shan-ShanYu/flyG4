#*差异表达分析--------------------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())#清空
Num = "012.2."
#### 1.剔除后的数据 ####
counts <- fread("/home/yuss/flyG4/result/qians.RNAseq/012.1.counts.txt") %>% as.data.frame()
rownames(counts) <- counts$geneid
counts <- counts[,-c(13,11)]
sample_name <- factor(colnames(counts))
metadata <- data.frame(sample_name)
metadata$treat <- as.factor(rep(c("DMSO", "TMP100", "TMP25","TMP50"), c(3, 3, 3,2)))

#### 2.DESeq2 TMPyP4 vs DMSO 火山图#### 
##1.构建矩阵
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ treat)
##2.dds标准化
dds <- DESeq(dds)
##3.获取标准化后的数据
normalized_counts <- counts(dds,normalized=T)
##4.有三个重复取平均
normalized_mean_counts = t(apply(normalized_counts, 1, function(a){tapply(a, metadata$treat, mean)}))

TMP100 <- na.omit(as.data.frame(results(dds, contrast = c("treat", "TMP100","DMSO"), cooksCutoff = FALSE))) #cooksCutoff = FALSE参数是指离群的基因不会用NA表示，之前没设置参数时，将会把离群值设为NA
TMP100$con <- normalized_mean_counts[match(rownames(TMP100),rownames(normalized_mean_counts)), 1]
TMP100$TMP100 <- normalized_mean_counts[match(rownames(TMP100),rownames(normalized_mean_counts)), 2]
TMP100$group <- ifelse(TMP100$pvalue<0.05&abs(TMP100$log2FoldChange)>=0.5,ifelse(TMP100$log2FoldChange>0.5,"Up","Down"),"No-sig")
table(TMP100$group)
# Down No-sig     Up 
# 609   6649    223
TMP100$geneid <- rownames(TMP100)
write.table(TMP100,file = paste0("/home/yuss/flyG4/result/qians.RNAseq/",Num,"DE.TMP100.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)

ggplot(TMP100,aes(x=log2FoldChange,y=-log10(pvalue))) +#黑白主题
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
  annotate(geom="text", x=-2.5, y=20, size=4,label=paste0("Down\n(N=",sum(TMP100$group=="Down"),")"))+
  annotate(geom="text", x=2, y=20, size=4,label=paste0("Up\n(N=",sum(TMP100$group=="Up"),")"))+
  guides(color="none") + #图例颜色
  ## 修改坐标轴
  xlab(bquote(Log[2]~Fold~Change))+
  ylab(bquote(Log[10]~P~Value(TMPyP4/DMSO))) +
  coord_cartesian(ylim = c(0,40),xlim = c(-4,3)) + 
  labs(title="TMPyP4 100umol/L") +
  annotate(geom = "text", x=2.5, y=38,size=4,
           label=paste0("Total:",nrow(TMP100))) 
# ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"DEKc.Phenvscon.VolcanoPlot.pdf"),
       device = "pdf",width = 4,height = 3)  

TMP25 <- na.omit(as.data.frame(results(dds, contrast = c("treat", "TMP25","DMSO"), cooksCutoff = FALSE))) #cooksCutoff = FALSE参数是指离群的基因不会用NA表示，之前没设置参数时，将会把离群值设为NA
TMP25$con <- normalized_mean_counts[match(rownames(TMP25),rownames(normalized_mean_counts)), 1]
TMP25$TMP25 <- normalized_mean_counts[match(rownames(TMP25),rownames(normalized_mean_counts)), 3]
TMP25$group <- ifelse(TMP25$pvalue<0.05&abs(TMP25$log2FoldChange)>=0.5,ifelse(TMP25$log2FoldChange>0.5,"Up","Down"),"No-sig")
table(TMP25$group)
# Down No-sig     Up 
# 349   6735    168
TMP25$geneid <- rownames(TMP25)
write.table(TMP25,file = paste0("/home/yuss/flyG4/result/qians.RNAseq/",Num,"DE.TMP25.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)
ggplot(TMP25,aes(x=log2FoldChange,y=-log10(pvalue))) +#黑白主题
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
  annotate(geom="text", x=-2.5, y=20, size=4,label=paste0("Down\n(N=",sum(TMP25$group=="Down"),")"))+
  annotate(geom="text", x=2, y=20, size=4,label=paste0("Up\n(N=",sum(TMP25$group=="Up"),")"))+
  guides(color="none") + #图例颜色
  ## 修改坐标轴
  xlab(bquote(Log[2]~Fold~Change))+
  ylab(bquote(Log[10]~P~Value(TMPyP4/DMSO))) +
  coord_cartesian(ylim = c(0,40),xlim = c(-4,3)) + 
  labs(title="TMPyP4 25umol/L") +
  annotate(geom = "text", x=2.5, y=38,size=4,
           label=paste0("Total:",nrow(TMP25))) 
# ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"DEKc.Phenvscon.VolcanoPlot.pdf"),
device = "pdf",width = 4,height = 3)  

TMP50 <- na.omit(as.data.frame(results(dds, contrast = c("treat", "TMP50","DMSO"), cooksCutoff = FALSE))) #cooksCutoff = FALSE参数是指离群的基因不会用NA表示，之前没设置参数时，将会把离群值设为NA
TMP50$con <- normalized_mean_counts[match(rownames(TMP50),rownames(normalized_mean_counts)), 1]
TMP50$TMP50 <- normalized_mean_counts[match(rownames(TMP50),rownames(normalized_mean_counts)), 4]
TMP50$group <- ifelse(TMP50$pvalue<0.05&abs(TMP50$log2FoldChange)>=0.5,ifelse(TMP50$log2FoldChange>0.5,"Up","Down"),"No-sig")
table(TMP50$group)
# Down No-sig     Up 
# 590   6561    101 
TMP50$geneid <- rownames(TMP50)
write.table(TMP50,file = paste0("/home/yuss/flyG4/result/qians.RNAseq/",Num,"DE.TMP50.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)
ggplot(TMP50,aes(x=log2FoldChange,y=-log10(pvalue))) +#黑白主题
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
  annotate(geom="text", x=-2.5, y=20, size=4,label=paste0("Down\n(N=",sum(TMP50$group=="Down"),")"))+
  annotate(geom="text", x=2, y=20, size=4,label=paste0("Up\n(N=",sum(TMP50$group=="Up"),")"))+
  guides(color="none") + #图例颜色
  ## 修改坐标轴
  xlab(bquote(Log[2]~Fold~Change))+
  ylab(bquote(Log[10]~P~Value(TMPyP4/DMSO))) +
  coord_cartesian(ylim = c(0,40),xlim = c(-4,3)) + 
  labs(title="TMPyP4 50umol/L") +
  annotate(geom = "text", x=2.5, y=38,size=4,
           label=paste0("Total:",nrow(TMP50))) 
# ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"DEKc.Phenvscon.VolcanoPlot.pdf"),
device = "pdf",width = 4,height = 3)  

#*TMPyP4_100----------------------------------------------------
#### 1.fisher含有eG4的基因在DEG还是非DEG上富集 #### 
rm(list = ls());gc();rm(list = ls())#清空
Num = "012.2."
T100 <- fread("/home/yuss/flyG4/result/qians.RNAseq/012.2.DE.TMP100.txt") %>% as.data.frame()
gene.s2 <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.gene.s2.bed") %>% as.data.frame()
T100$num <- gene.s2[match(T100$geneid,gene.s2$V4),7]
T100$s2.eG4 <- ifelse(T100$num==0,"no eG4","eG4")
table(T100$s2.eG4,T100$group)
# Down No-sig   Up
# eG4     352   1352   27
# no eG4  257   5297  196
# DEG eG4= 352+27=379
# DEG no eG4=257+196=453
data <- matrix(round((c(379,1352,453,5297)),0),nrow = 2)
fisher.test(data) # p-value < 2.2e-16
data <- matrix(round((c(352,1352,257,5297)),0),nrow = 2)
fisher.test(data) # p-value < 2.2e-16
data <- matrix(round((c(27,1352,196,5297)),0),nrow = 2)
fisher.test(data) # p-value = 0.002106 #上调基因缺失eG4

#### 2.TMPyP4_100基因上eG4的信号(密度) ####
T100$chr <- gene.s2[match(T100$geneid,gene.s2$V4),1]
T100$start <- gene.s2[match(T100$geneid,gene.s2$V4),2]
T100$end <- gene.s2[match(T100$geneid,gene.s2$V4),3]
T100$length <- T100$end-T100$start
T100$density <- T100$num/T100$length

T100$group <- factor(T100$group,levels = c("Down","No-sig","Up"))
my_comparisons = list(c("No-sig","Down"),c("No-sig","Up"),c("Down","Up"))

ggplot(data = T100,aes(x=group,y=density*1000,fill=group)) +
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
tapply(T100$density, T100$group, median)

T100$numtype <- ifelse(T100$num > 0,"have","no")
a <- as.data.frame(table(T100$group, T100$numtype))
wide_a <- spread(a, Var2, Freq)
wide_a$sum <- wide_a$have + wide_a$no
wide_a$have.ratio <- wide_a$have/wide_a$sum
wide_a$no.ratio <- wide_a$no/wide_a$sum
wide_a$Var1 <- factor(wide_a$Var1,levels = c("Down","No-sig","Up"))

library(ggsci)
ggplot(wide_a, aes(x=Var1,y=have.ratio*100,fill=Var1)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(0.7),width = 0.7,color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值
  coord_cartesian(ylim = c(0,80)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(name = '',
                    values = c('Up'='#D6604D','No-sig'='grey','Down'='#74ADD1'), #手动设置颜色时调整颜色的因子顺序
  )+
  # geom_text(aes(label=have.ratio*100,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题(左下边框，没有网格线)
  ylab("(%) Gene with eG4") + ##通过bquote函数给图标签添加上下标
  geom_signif(y_position=62, xmin=1, xmax=2,
              annotation=c("2.2e-16"),tip_length=0)+
  geom_signif(y_position=25, xmin=2, xmax=3,
              annotation=c("0.002"),tip_length=0)+
  geom_signif(y_position=68, xmin=1, xmax=3,
              annotation=c("2.2e-16"),tip_length=0)+
  theme(plot.title = element_text(size = 14, face = "plain"),
        axis.title.y = element_text(size = 14), ##y坐标轴标题字体大小
        axis.title.x = element_blank(), ##删除x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.position = "none") +
  labs(title="TMPyP4 100umol/L")
# ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"KcPhenDEG.PercentagewitheG4.pdf"),
#        device = "pdf",width = 3.2,height = 3.4)
data <- matrix((c(352,1352,257,5297)),nrow = 2)
fisher.test(data) #p-value < 2.2e-16

data <- matrix((c(27,1352,196,5297)),nrow = 2)
fisher.test(data) #p-value = 0.002106

data <- matrix((c(352,27,257,196)),nrow = 2)
fisher.test(data) #p-value < 2.2e-16

#### 3.TMPyP4_100基因的Promoter区eG4的信号 ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "012.2."
promoter.s2G4 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter2000.s2.bed") %>% as.data.frame()
t100 <- fread("/home/yuss/flyG4/result/qians.RNAseq/012.2.DE.TMP100.txt") %>% as.data.frame()
table(t100$group)
t100$promoter.num <- promoter.s2G4[match(t100$geneid,promoter.s2G4$V4),6]
t100$group <- factor(t100$group,levels = c("Down","No-sig","Up"))
my_comparisons = list(c("No-sig","Down"),c("No-sig","Up"))

ggplot(data = t100,aes(x=group,y=promoter.num,fill=group)) +
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
  ylab("Number of Promoter eG4") +
  coord_cartesian(ylim = c(0,0.3)) +
  labs(title="TMPyp4")

t100$numtype <- ifelse(t100$promoter.num > 0,"have","no")
a <- as.data.frame(table(t100$group, t100$numtype))
wide_a <- spread(a, Var2, Freq)
wide_a$sum <- wide_a$have + wide_a$no
wide_a$have.ratio <- wide_a$have/wide_a$sum
wide_a$no.ratio <- wide_a$no/wide_a$sum
wide_a$Var1 <- factor(wide_a$Var1,levels = c("Down","No-sig","Up"))

library(ggsci)
ggplot(wide_a, aes(x=Var1,y=have.ratio*100,fill=Var1)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(0.7),width = 0.7,color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值
  coord_cartesian(ylim = c(0,40)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(name = '',
                    values = c('Up'='#D6604D','No-sig'='grey','Down'='#74ADD1'), #手动设置颜色时调整颜色的因子顺序
  )+
  # geom_text(aes(label=have.ratio*100,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题(左下边框，没有网格线)
  ylab("(%) Gene with Promoter eG4") + ##通过bquote函数给图标签添加上下标
  geom_signif(y_position=33, xmin=1, xmax=2,
              annotation=c("1.69e-07"),tip_length=0)+
  geom_signif(y_position=24, xmin=2, xmax=3,
              annotation=c("0.098"),tip_length=0)+
  geom_signif(y_position=37, xmin=1, xmax=3,
              annotation=c("2.47e-05"),tip_length=0)+
  theme(plot.title = element_text(size = 14, face = "plain"),
        axis.title.y = element_text(size = 14), ##y坐标轴标题字体大小
        axis.title.x = element_blank(), ##删除x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.position = "none") +
  labs(title="TMPyP4 100umol/L")
# ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"KcPhenDEG.PercentagewitheG4.pdf"),
#        device = "pdf",width = 3.2,height = 3.4)
data <- matrix((c(192,1454,417,5195)),nrow = 2)
fisher.test(data) #p-value = 1.688e-07

data <- matrix((c(38,1454,185,5195)),nrow = 2)
fisher.test(data) #p-value = 0.0982

data <- matrix((c(38,192,185,417)),nrow = 2)
fisher.test(data) #p-value = 2.473e-05

#*TMPyP4_25----------------------------------------------------
#### 1.fisher含有eG4的基因在DEG还是非DEG上富集 #### 
rm(list = ls());gc();rm(list = ls())#清空
Num = "012.2."
t25 <- fread("/home/yuss/flyG4/result/qians.RNAseq/012.2.DE.TMP25.txt") %>% as.data.frame()
gene.s2 <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.gene.s2.bed") %>% as.data.frame()
t25$num <- gene.s2[match(t25$geneid,gene.s2$V4),7]
t25$s2.eG4 <- ifelse(t25$num==0,"no eG4","eG4")
table(t25$s2.eG4,t25$group)
# Down No-sig   Up
# eG4     203   1437   30
# no eG4  146   5298  138
# DEG eG4= 203+30=233
# DEG no eG4=146+138=284
data <- matrix(round((c(233,1437,284,5298)),0),nrow = 2)
fisher.test(data) # p-value < 2.2e-16
data <- matrix(round((c(203,1437,146,5298)),0),nrow = 2)
fisher.test(data) # p-value < 2.2e-16
data <- matrix(round((c(30,1437,138,5298)),0),nrow = 2)
fisher.test(data) # p-value = 0.2951 #上调基因不富集也不缺失eG4

#### 2.TMPyP4_25基因上eG4的信号(密度) ####
t25$chr <- gene.s2[match(t25$geneid,gene.s2$V4),1]
t25$start <- gene.s2[match(t25$geneid,gene.s2$V4),2]
t25$end <- gene.s2[match(t25$geneid,gene.s2$V4),3]
t25$length <- t25$end-t25$start
t25$density <- t25$num/t25$length

t25$group <- factor(t25$group,levels = c("Down","No-sig","Up"))
my_comparisons = list(c("No-sig","Down"),c("No-sig","Up"),c("Down","Up"))

ggplot(data = t25,aes(x=group,y=density*1000,fill=group)) +
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
tapply(t25$density, t25$group, median)

t25$numtype <- ifelse(t25$num > 0,"have","no")
a <- as.data.frame(table(t25$group, t25$numtype))
wide_a <- spread(a, Var2, Freq)
wide_a$sum <- wide_a$have + wide_a$no
wide_a$have.ratio <- wide_a$have/wide_a$sum
wide_a$no.ratio <- wide_a$no/wide_a$sum
wide_a$Var1 <- factor(wide_a$Var1,levels = c("Down","No-sig","Up"))

library(ggsci)
ggplot(wide_a, aes(x=Var1,y=have.ratio*100,fill=Var1)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(0.7),width = 0.7,color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值
  coord_cartesian(ylim = c(0,80)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(name = '',
                    values = c('Up'='#D6604D','No-sig'='grey','Down'='#74ADD1'), #手动设置颜色时调整颜色的因子顺序
  )+
  # geom_text(aes(label=have.ratio*100,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题(左下边框，没有网格线)
  ylab("(%) Gene with eG4") + ##通过bquote函数给图标签添加上下标
  geom_signif(y_position=62, xmin=1, xmax=2,
              annotation=c("2.2e-16"),tip_length=0)+
  geom_signif(y_position=25, xmin=2, xmax=3,
              annotation=c("0.30"),tip_length=0)+
  geom_signif(y_position=68, xmin=1, xmax=3,
              annotation=c("2.2e-16"),tip_length=0)+
  theme(plot.title = element_text(size = 14, face = "plain"),
        axis.title.y = element_text(size = 14), ##y坐标轴标题字体大小
        axis.title.x = element_blank(), ##删除x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.position = "none") +
  labs(title="TMPyP4 25umol/L")
# ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"KcPhenDEG.PercentagewitheG4.pdf"),
#        device = "pdf",width = 3.2,height = 3.4)
data <- matrix((c(203,1437,146,5298)),nrow = 2)
fisher.test(data) #p-value < 2.2e-16

data <- matrix((c(30,1437,138,5298)),nrow = 2)
fisher.test(data) #p-value = 0.2951

data <- matrix((c(203,30,146,138)),nrow = 2)
fisher.test(data) #p-value < 2.2e-16

#### 3.TMPyP4_25基因的Promoter区eG4的信号 ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "012.2."
promoter.s2G4 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter2000.s2.bed") %>% as.data.frame()
t25 <- fread("/home/yuss/flyG4/result/qians.RNAseq/012.2.DE.TMP25.txt") %>% as.data.frame()
table(t25$group)
t25$promoter.num <- promoter.s2G4[match(t25$geneid,promoter.s2G4$V4),6]
t25$group <- factor(t25$group,levels = c("Down","No-sig","Up"))
my_comparisons = list(c("No-sig","Down"),c("No-sig","Up"))

ggplot(data = t25,aes(x=group,y=promoter.num,fill=group)) +
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
  ylab("Number of Promoter eG4") +
  coord_cartesian(ylim = c(0,0.3)) +
  labs(title="TMPyp4")

t25$numtype <- ifelse(t25$promoter.num > 0,"have","no")
a <- as.data.frame(table(t25$group, t25$numtype))
wide_a <- spread(a, Var2, Freq)
wide_a$sum <- wide_a$have + wide_a$no
wide_a$have.ratio <- wide_a$have/wide_a$sum
wide_a$no.ratio <- wide_a$no/wide_a$sum
wide_a$Var1 <- factor(wide_a$Var1,levels = c("Down","No-sig","Up"))

library(ggsci)
ggplot(wide_a, aes(x=Var1,y=have.ratio*100,fill=Var1)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(0.7),width = 0.7,color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值
  coord_cartesian(ylim = c(0,40)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(name = '',
                    values = c('Up'='#D6604D','No-sig'='grey','Down'='#74ADD1'), #手动设置颜色时调整颜色的因子顺序
  )+
  # geom_text(aes(label=have.ratio*100,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题(左下边框，没有网格线)
  ylab("(%) Gene with Promoter eG4") + ##通过bquote函数给图标签添加上下标
  geom_signif(y_position=34, xmin=1, xmax=2,
              annotation=c("3.41e-06"),tip_length=0)+
  geom_signif(y_position=25, xmin=2, xmax=3,
              annotation=c("0.51"),tip_length=0)+
  geom_signif(y_position=37, xmin=1, xmax=3,
              annotation=c("0.040"),tip_length=0)+
  theme(plot.title = element_text(size = 14, face = "plain"),
        axis.title.y = element_text(size = 14), ##y坐标轴标题字体大小
        axis.title.x = element_blank(), ##删除x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.position = "none") +
  labs(title="TMPyP4 25umol/L")
# ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"KcPhenDEG.PercentagewitheG4.pdf"),
#        device = "pdf",width = 3.2,height = 3.4)
data <- matrix((c(115,1467,234,5268)),nrow = 2)
fisher.test(data) #p-value = 3.409e-06

data <- matrix((c(40,1467,128,5268)),nrow = 2)
fisher.test(data) #p-value = 0.5094

data <- matrix((c(115,40,234,128)),nrow = 2)
fisher.test(data) #p-value = 0.04019

#*TMPyP4_50----------------------------------------------------
#### 1.fisher含有eG4的基因在DEG还是非DEG上富集 #### 
rm(list = ls());gc();rm(list = ls())#清空
Num = "012.2."
t50 <- fread("/home/yuss/flyG4/result/qians.RNAseq/012.2.DE.TMP50.txt") %>% as.data.frame()
gene.s2 <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.gene.s2.bed") %>% as.data.frame()
t50$num <- gene.s2[match(t50$geneid,gene.s2$V4),7]
t50$s2.eG4 <- ifelse(t50$num==0,"no eG4","eG4")
table(t50$s2.eG4,t50$group)
# Down No-sig   Up
# eG4     334   1323   13
# no eG4  256   5238   88
# DEG eG4= 334+13=347
# DEG no eG4=256+88=344
data <- matrix(round((c(347,1323,344,5238)),0),nrow = 2)
fisher.test(data) # p-value < 2.2e-16
data <- matrix(round((c(334,1323,256,5238)),0),nrow = 2)
fisher.test(data) # p-value < 2.2e-16
data <- matrix(round((c(13,1323,88,5238)),0),nrow = 2)
fisher.test(data) # p-value = 0.07882 #上调基因

#### 2.TMPyP4_50基因上eG4的信号(密度) ####
t50$chr <- gene.s2[match(t50$geneid,gene.s2$V4),1]
t50$start <- gene.s2[match(t50$geneid,gene.s2$V4),2]
t50$end <- gene.s2[match(t50$geneid,gene.s2$V4),3]
t50$length <- t50$end-t50$start
t50$density <- t50$num/t50$length

t50$group <- factor(t50$group,levels = c("Down","No-sig","Up"))
my_comparisons = list(c("No-sig","Down"),c("No-sig","Up"),c("Down","Up"))

ggplot(data = t50,aes(x=group,y=density*1000,fill=group)) +
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
tapply(t50$density, t50$group, median)

t50$numtype <- ifelse(t50$num > 0,"have","no")
a <- as.data.frame(table(t50$group, t50$numtype))
wide_a <- spread(a, Var2, Freq)
wide_a$sum <- wide_a$have + wide_a$no
wide_a$have.ratio <- wide_a$have/wide_a$sum
wide_a$no.ratio <- wide_a$no/wide_a$sum
wide_a$Var1 <- factor(wide_a$Var1,levels = c("Down","No-sig","Up"))

library(ggsci)
ggplot(wide_a, aes(x=Var1,y=have.ratio*100,fill=Var1)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(0.7),width = 0.7,color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值
  coord_cartesian(ylim = c(0,80)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(name = '',
                    values = c('Up'='#D6604D','No-sig'='grey','Down'='#74ADD1'), #手动设置颜色时调整颜色的因子顺序
  )+
  # geom_text(aes(label=have.ratio*100,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题(左下边框，没有网格线)
  ylab("(%) Gene with eG4") + ##通过bquote函数给图标签添加上下标
  geom_signif(y_position=61, xmin=1, xmax=2,
              annotation=c("2.2e-16"),tip_length=0)+
  geom_signif(y_position=25, xmin=2, xmax=3,
              annotation=c("0.079"),tip_length=0)+
  geom_signif(y_position=68, xmin=1, xmax=3,
              annotation=c("2.2e-16"),tip_length=0)+
  theme(plot.title = element_text(size = 14, face = "plain"),
        axis.title.y = element_text(size = 14), ##y坐标轴标题字体大小
        axis.title.x = element_blank(), ##删除x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.position = "none") +
  labs(title="TMPyP4 50umol/L")
# ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"KcPhenDEG.PercentagewitheG4.pdf"),
#        device = "pdf",width = 3.2,height = 3.4)
data <- matrix((c(334,1323,256,5238)),nrow = 2)
fisher.test(data) #p-value < 2.2e-16

data <- matrix((c(13,1323,88,5238)),nrow = 2)
fisher.test(data) #p-value = 0.07882

data <- matrix((c(334,13,256,88)),nrow = 2)
fisher.test(data) #p-value < 2.2e-16

#### 3.TMPyP4_50基因的Promoter区eG4的信号 ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "012.2."
promoter.s2G4 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter2000.s2.bed") %>% as.data.frame()
t50 <- fread("/home/yuss/flyG4/result/qians.RNAseq/012.2.DE.TMP50.txt") %>% as.data.frame()
table(t50$group)
t50$promoter.num <- promoter.s2G4[match(t50$geneid,promoter.s2G4$V4),6]
t50$group <- factor(t50$group,levels = c("Down","No-sig","Up"))
my_comparisons = list(c("No-sig","Down"),c("No-sig","Up"))

ggplot(data = t50,aes(x=group,y=promoter.num,fill=group)) +
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
  ylab("Number of Promoter eG4") +
  coord_cartesian(ylim = c(0,0.3)) +
  labs(title="TMPyp4")

t50$numtype <- ifelse(t50$promoter.num > 0,"have","no")
a <- as.data.frame(table(t50$group, t50$numtype))
wide_a <- spread(a, Var2, Freq)
wide_a$sum <- wide_a$have + wide_a$no
wide_a$have.ratio <- wide_a$have/wide_a$sum
wide_a$no.ratio <- wide_a$no/wide_a$sum
wide_a$Var1 <- factor(wide_a$Var1,levels = c("Down","No-sig","Up"))

library(ggsci)
ggplot(wide_a, aes(x=Var1,y=have.ratio*100,fill=Var1)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(0.7),width = 0.7,color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值
  coord_cartesian(ylim = c(0,40)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(name = '',
                    values = c('Up'='#D6604D','No-sig'='grey','Down'='#74ADD1'), #手动设置颜色时调整颜色的因子顺序
  )+
  # geom_text(aes(label=have.ratio*100,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题(左下边框，没有网格线)
  ylab("(%) Gene with Promoter eG4") + ##通过bquote函数给图标签添加上下标
  geom_signif(y_position=33, xmin=1, xmax=2,
              annotation=c("1.37e-09"),tip_length=0)+
  geom_signif(y_position=24, xmin=2, xmax=3,
              annotation=c("0.067"),tip_length=0)+
  geom_signif(y_position=37, xmin=1, xmax=3,
              annotation=c("5.90e-05"),tip_length=0)+
  theme(plot.title = element_text(size = 14, face = "plain"),
        axis.title.y = element_text(size = 14), ##y坐标轴标题字体大小
        axis.title.x = element_blank(), ##删除x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.position = "none") +
  labs(title="TMPyP4 50umol/L")
# ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"KcPhenDEG.PercentagewitheG4.pdf"),
#        device = "pdf",width = 3.2,height = 3.4)
data <- matrix((c(194,1414,396,5147)),nrow = 2)
fisher.test(data) #p-value =1.371e-09

data <- matrix((c(14,1414,87,5147)),nrow = 2)
fisher.test(data) #p-value = 0.06651

data <- matrix((c(194,14,396,87)),nrow = 2)
fisher.test(data) #p-value = 5.904e-05

#*含有eG4和不含有eG4的基因的差异倍数 log2FC---------------------------------------------------------------------------
#### 1.TMPyP4_100 基因含有eG4和no eG4差异倍数  ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "012.2."
t100 <- fread("/home/yuss/flyG4/result/qians.RNAseq/012.2.DE.TMP100.txt") %>% as.data.frame()
gene.s2 <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.gene.s2.bed") %>% as.data.frame()
t100$num <- gene.s2[match(t100$geneid,gene.s2$V4),7]
t100$s2.eG4 <- ifelse(t100$num==0,"no eG4","eG4")
t100$s2.eG4 <- factor(t100$s2.eG4,levels = c("no eG4","eG4"))
my_comparisons = list(c("eG4","no eG4"))
ggplot(data = t100,aes(x=s2.eG4,y=log2FoldChange,fill=s2.eG4,color=s2.eG4)) +
  geom_violin(position = position_dodge(width = 1), scale = 'width') +
  geom_boxplot(position = position_dodge(width = 1), outlier.size = 0.6, width = 0.2, show.legend = FALSE) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(3.5,5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(-5,5)) +
  scale_fill_manual(values = c("#818181","#D08725")) +
  scale_color_manual(values = c("#5b5b5b","#8f5d19")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2](TMPyP4/DMSO))) +
  labs(title="TMPyP4 100umol/L") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") 
# ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"S2PhenlogFC.witheG4noeG4.pdf"),
       device = "pdf",width = 2.7,height = 2.7)

tapply(t100$log2FoldChange,t100$s2.eG4, mean)
tapply(t100$log2FoldChange,t100$s2.eG4, median)

sample_data <- t100[t100$s2.eG4=="eG4",2]
# 单样本Wilcoxon检验(看小于0是否有显著性)
wilcox.test(sample_data, mu = 0,alternative = "less") #p-value < 2.2e-16
sample_data <- t100[t100$s2.eG4=="no eG4",2]
# 单样本Wilcoxon检验(看大于0是否有显著性)
wilcox.test(sample_data, mu = 0,alternative = "greater") #p-value = 2.079e-12 显著大于0

#### 2.TMPyP4_100 启动子含有eG4和no eG4差异倍数  ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "012.2."
t100 <- fread("/home/yuss/flyG4/result/qians.RNAseq/012.2.DE.TMP100.txt") %>% as.data.frame()
promoter.s2G4 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter2000.s2.bed") %>% as.data.frame()
t100$pro.num <- promoter.s2G4[match(t100$geneid,promoter.s2G4$V4),6]
t100$s2.eG4 <- ifelse(t100$pro.num==0,"no eG4","eG4")
t100$s2.eG4 <- factor(t100$s2.eG4,levels = c("no eG4","eG4"))
my_comparisons = list(c("eG4","no eG4"))
ggplot(data = t100,aes(x=s2.eG4,y=log2FoldChange,fill=s2.eG4,color=s2.eG4)) +
  geom_violin(position = position_dodge(width = 1), scale = 'width') +
  geom_boxplot(position = position_dodge(width = 1), outlier.size = 0.6, width = 0.2, show.legend = FALSE) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(3.5,5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(-5,5)) +
  scale_fill_manual(values = c("#818181","#D08725")) +
  scale_color_manual(values = c("#5b5b5b","#8f5d19")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2](TMPyP4/DMSO))) +
  labs(title="TMPyP4 100umol/L (Promter eG4)") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") 
# ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"S2PhenlogFC.withPromotereG4noeG4.pdf"),
       device = "pdf",width = 2.7,height = 2.7)

tapply(t100$log2FoldChange,t100$s2.eG4, mean)
tapply(t100$log2FoldChange,t100$s2.eG4, median)

sample_data <- t100[t100$s2.eG4=="eG4",2]
# 单样本Wilcoxon检验(看大于0是否有显著性)
wilcox.test(sample_data, mu = 0,alternative = "less") #p-value = 1.332e-11
sample_data <- t100[t100$s2.eG4=="no eG4",2]
# 单样本Wilcoxon检验(看小于0是否有显著性)
wilcox.test(sample_data, mu = 0,alternative = "less") #p-value = 0.03727

#### 3.TMPyP4_25 基因含有eG4和no eG4差异倍数  ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "012.2."
t25 <- fread("/home/yuss/flyG4/result/qians.RNAseq/012.2.DE.TMP25.txt") %>% as.data.frame()
gene.s2 <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.gene.s2.bed") %>% as.data.frame()
t25$num <- gene.s2[match(t25$geneid,gene.s2$V4),7]
t25$s2.eG4 <- ifelse(t25$num==0,"no eG4","eG4")
t25$s2.eG4 <- factor(t25$s2.eG4,levels = c("no eG4","eG4"))
my_comparisons = list(c("eG4","no eG4"))
ggplot(data = t25,aes(x=s2.eG4,y=log2FoldChange,fill=s2.eG4,color=s2.eG4)) +
  geom_violin(position = position_dodge(width = 1), scale = 'width') +
  geom_boxplot(position = position_dodge(width = 1), outlier.size = 0.6, width = 0.2, show.legend = FALSE) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(3.5,5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(-5,5)) +
  scale_fill_manual(values = c("#818181","#D08725")) +
  scale_color_manual(values = c("#5b5b5b","#8f5d19")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2](TMPyP4/DMSO))) +
  labs(title="TMPyP4 25umol/L") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") 
# ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"S2PhenlogFC.witheG4noeG4.pdf"),
device = "pdf",width = 2.7,height = 2.7)

tapply(t25$log2FoldChange,t25$s2.eG4, mean)
tapply(t25$log2FoldChange,t25$s2.eG4, median)

sample_data <- t25[t25$s2.eG4=="eG4",2]
# 单样本Wilcoxon检验(看小于0是否有显著性)
wilcox.test(sample_data, mu = 0,alternative = "less") #p-value < 2.2e-16
sample_data <- t25[t25$s2.eG4=="no eG4",2]
# 单样本Wilcoxon检验(看大于0是否有显著性)
wilcox.test(sample_data, mu = 0,alternative = "greater") #p-value = 1.528e-06 显著大于0

#### 4.TMPyP4_25 启动子含有eG4和no eG4差异倍数  ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "012.2."
t25 <- fread("/home/yuss/flyG4/result/qians.RNAseq/012.2.DE.TMP25.txt") %>% as.data.frame()
promoter.s2G4 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter2000.s2.bed") %>% as.data.frame()
t25$pro.num <- promoter.s2G4[match(t25$geneid,promoter.s2G4$V4),6]
t25$s2.eG4 <- ifelse(t25$pro.num==0,"no eG4","eG4")
t25$s2.eG4 <- factor(t25$s2.eG4,levels = c("no eG4","eG4"))
my_comparisons = list(c("eG4","no eG4"))
ggplot(data = t25,aes(x=s2.eG4,y=log2FoldChange,fill=s2.eG4,color=s2.eG4)) +
  geom_violin(position = position_dodge(width = 1), scale = 'width') +
  geom_boxplot(position = position_dodge(width = 1), outlier.size = 0.6, width = 0.2, show.legend = FALSE) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(3.5,5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(-5,5)) +
  scale_fill_manual(values = c("#818181","#D08725")) +
  scale_color_manual(values = c("#5b5b5b","#8f5d19")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2](TMPyP4/DMSO))) +
  labs(title="TMPyP4 25umol/L (Promter eG4)") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") 
# ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"S2PhenlogFC.withPromotereG4noeG4.pdf"),
device = "pdf",width = 2.7,height = 2.7)

tapply(t25$log2FoldChange,t25$s2.eG4, mean)
tapply(t25$log2FoldChange,t25$s2.eG4, median)

sample_data <- t25[t25$s2.eG4=="eG4",2]
# 单样本Wilcoxon检验(看大于0是否有显著性)
wilcox.test(sample_data, mu = 0,alternative = "less") #p-value = 2.304e-08
sample_data <- t25[t25$s2.eG4=="no eG4",2]
# 单样本Wilcoxon检验(看小于0是否有显著性)
wilcox.test(sample_data, mu = 0,alternative = "less") #p-value = 0.02589

#### 5.TMPyP4_50 基因含有eG4和no eG4差异倍数  ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "012.2."
t50 <- fread("/home/yuss/flyG4/result/qians.RNAseq/012.2.DE.TMP50.txt") %>% as.data.frame()
gene.s2 <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.gene.s2.bed") %>% as.data.frame()
t50$num <- gene.s2[match(t50$geneid,gene.s2$V4),7]
t50$s2.eG4 <- ifelse(t50$num==0,"no eG4","eG4")
t50$s2.eG4 <- factor(t50$s2.eG4,levels = c("no eG4","eG4"))
my_comparisons = list(c("eG4","no eG4"))
ggplot(data = t50,aes(x=s2.eG4,y=log2FoldChange,fill=s2.eG4,color=s2.eG4)) +
  geom_violin(position = position_dodge(width = 1), scale = 'width') +
  geom_boxplot(position = position_dodge(width = 1), outlier.size = 0.6, width = 0.2, show.legend = FALSE) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(3.5,5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(-5,5)) +
  scale_fill_manual(values = c("#818181","#D08725")) +
  scale_color_manual(values = c("#5b5b5b","#8f5d19")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2](TMPyP4/DMSO))) +
  labs(title="TMPyP4 50umol/L") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") 
# ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"S2PhenlogFC.witheG4noeG4.pdf"),
device = "pdf",width = 2.7,height = 2.7)

tapply(t50$log2FoldChange,t50$s2.eG4, mean)
tapply(t50$log2FoldChange,t50$s2.eG4, median)

sample_data <- t50[t50$s2.eG4=="eG4",2]
# 单样本Wilcoxon检验(看小于0是否有显著性)
wilcox.test(sample_data, mu = 0,alternative = "less") #p-value < 2.2e-16
sample_data <- t50[t50$s2.eG4=="no eG4",2]
# 单样本Wilcoxon检验(看大于0是否有显著性)
wilcox.test(sample_data, mu = 0,alternative = "greater") #p-value = 1.796e-07 显著大于0

#### 6.TMPyP4_50 启动子含有eG4和no eG4差异倍数  ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "012.2."
t50 <- fread("/home/yuss/flyG4/result/qians.RNAseq/012.2.DE.TMP50.txt") %>% as.data.frame()
promoter.s2G4 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter2000.s2.bed") %>% as.data.frame()
t50$pro.num <- promoter.s2G4[match(t50$geneid,promoter.s2G4$V4),6]
t50$s2.eG4 <- ifelse(t50$pro.num==0,"no eG4","eG4")
t50$s2.eG4 <- factor(t50$s2.eG4,levels = c("no eG4","eG4"))
my_comparisons = list(c("eG4","no eG4"))
ggplot(data = t50,aes(x=s2.eG4,y=log2FoldChange,fill=s2.eG4,color=s2.eG4)) +
  geom_violin(position = position_dodge(width = 1), scale = 'width') +
  geom_boxplot(position = position_dodge(width = 1), outlier.size = 0.6, width = 0.2, show.legend = FALSE) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(3.5,5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(-5,5)) +
  scale_fill_manual(values = c("#818181","#D08725")) +
  scale_color_manual(values = c("#5b5b5b","#8f5d19")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2](TMPyP4/DMSO))) +
  labs(title="TMPyP4 50umol/L (Promter eG4)") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") 
# ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"S2PhenlogFC.withPromotereG4noeG4.pdf"),
device = "pdf",width = 2.7,height = 2.7)

tapply(t50$log2FoldChange,t50$s2.eG4, mean)
tapply(t50$log2FoldChange,t50$s2.eG4, median)

sample_data <- t50[t50$s2.eG4=="eG4",2]
# 单样本Wilcoxon检验(看大于0是否有显著性)
wilcox.test(sample_data, mu = 0,alternative = "less") #p-value < 2.2e-16
sample_data <- t50[t50$s2.eG4=="no eG4",2]
# 单样本Wilcoxon检验(看小于0是否有显著性)
wilcox.test(sample_data, mu = 0,alternative = "less") #p-value = 5.61e-05
