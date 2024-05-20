rm(list = ls());gc();rm(list = ls())#清空
setwd("/home/yuss/flyG4/script")
Num = "001.5."
#### genome annotation (mergeG4和non_eG4)####
#*eG4基因组注释----------------------------------------------------------------------------------
##创建转录组数据库（TxDb)
##直接使用TxDb包中提供的函数来构建TxDb
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
Dmel.txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
##也可以用makeTxDbFromGFF()函数从GFF文件中创建TxDb   Dmel.txdb <-makeTxDbFromGFF("/home/qians/Quadruplex/Input/Ref/Fly/dmel-all-r6.19.gff",format = "gff3")
library(ChIPseeker)
merge.anno <- annotatePeak("/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.merge.chr.bed",
                           tssRegion = c(-2000,2000),TxDb = Dmel.txdb)
merge.anno
## 可视化
plotAnnoBar(merge.anno)
plotAnnoPie(merge.anno)
vennpie(merge.anno)

df = merge.anno@annoStat
class(df$Feature) ##factor
df$Feature %<>% as.character()
df[nrow(df)+1,] = c("Exon",df[df$Feature=="1st Exon",2]+df[df$Feature=="Other Exon",2])
df[nrow(df)+1,] = c("Intron", as.numeric(df[df$Feature=="1st Intron",2])+as.numeric(df[df$Feature=="Other Intron",2]))
df[nrow(df)+1,] = c("Promoter",as.numeric(df[df$Feature=="Promoter (<=1kb)",2])+as.numeric(df[df$Feature=="Promoter (1-2kb)",2]))
df %<>% dplyr::filter(Feature %in% c("Promoter","5' UTR","3' UTR","Exon","Intron","Distal Intergenic","Downstream (<=300)"))
df$Frequency %<>% as.numeric()
df
df %<>% arrange(-Frequency)
## df 中的 Feature 列中的特定字符串 "Downstream (<=300)" 替换为 "Downstream"，并将修改后的结果存回 df$Feature 列,fixed = TRUE 参数表示使用字符串匹配而不是正则表达式匹配,接着把这列设置为因子，rev(.) 表示将列的因子水平逆序排列。
df$Feature %<>% gsub("Downstream (<=300)","Downstream",.,fixed = TRUE) %>% factor(.,levels = rev(.))
ggplot(df,aes(Feature,Frequency))+
  geom_bar(stat = "identity",width = 0.75,fill="#F0988C") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 14),axis.text.x = element_text(size = 12,colour = "black"),
    axis.title.y = element_text(size = 14),axis.text.y = element_text(size = 12,colour = "black"),
    plot.margin = margin(t = 10,  # 顶部边缘距离
                         r = 10,  # 右边边缘距离
                         b = 10,  # 底部边缘距离
                         l = 10)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,50),position = "right") + ##y轴坐标轴在右边
  ylab("Frequancy (%)") +
  xlab("") +
  annotate(x=7.6, xend = 7.6, y=0, yend = 50, colour="black", lwd=0.75,geom="segment") +
  coord_flip()
  
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"LocationMergeG4.Dmel.pdf"),
        device = "pdf",width = 3,height = 5)
#* non_eG4基因组注释------------------------------------------------------------------------------------                  
non_eG4.anno <- annotatePeak("/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.non_eG4.chr.bed",
                           tssRegion = c(-2000,2000),TxDb = Dmel.txdb)
non_eG4.anno
##可视化
plotAnnoBar(non_eG4.anno)
plotAnnoPie(non_eG4.anno)
vennpie(non_eG4.anno)

df = non_eG4.anno@annoStat
class(df$Feature) ##factor
df$Feature %<>% as.character()
df[nrow(df)+1,] = c("Exon",df[df$Feature=="1st Exon",2]+df[df$Feature=="Other Exon",2])
df[nrow(df)+1,] = c("Intron", as.numeric(df[df$Feature=="1st Intron",2])+as.numeric(df[df$Feature=="Other Intron",2]))
df[nrow(df)+1,] = c("Promoter",as.numeric(df[df$Feature=="Promoter (<=1kb)",2])+as.numeric(df[df$Feature=="Promoter (1-2kb)",2]))
df %<>% dplyr::filter(Feature %in% c("Promoter","5' UTR","3' UTR","Exon","Intron","Distal Intergenic","Downstream (<=300)"))
df$Frequency %<>% as.numeric()
df
## df 中的 Feature 列中的特定字符串 "Downstream (<=300)" 替换为 "Downstream"，并将修改后的结果存回 df$Feature 列,fixed = TRUE 参数表示使用字符串匹配而不是正则表达式匹配,接着把这列设置为因子，rev(.) 表示将列的因子水平逆序排列。
df$Feature %<>% gsub("Downstream (<=300)","Downstream",.,fixed = TRUE) 
df$Feature <- factor(df$Feature,levels = c("Downstream","5' UTR","3' UTR","Exon","Distal Intergenic","Intron","Promoter"))
#df %<>% arrange(-Frequency)

ggplot(df,aes(Feature,Frequency))+
  geom_bar(stat = "identity",width = 0.75,fill="#A1A9D0") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 14),axis.text.x = element_text(size = 12,colour = "black"),
    axis.title.y = element_text(size = 14),axis.text.y = element_text(size = 12,colour = "black"),
    plot.margin = margin(t = 10,  # 顶部边缘距离
                         r = 10,  # 右边边缘距离
                         b = 10,  # 底部边缘距离
                         l = 10)) +
  scale_y_continuous(expand = c(0,0),limits = c(0,50),position = "right") + ##y轴坐标轴在右边
  ylab("Frequancy (%)") +
  xlab("") +
  annotate(x=7.6, xend = 7.6, y=0, yend = 50, colour="black", lwd=0.75,geom="segment") +
  coord_flip()
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"Locationnon_eG4.G4.Dmel.pdf"),
       device = "pdf",width = 3,height = 5)

rm(list = ls());gc();rm(list = ls())#清空
Num = "001.5."
#### Length distribution ####
mergeG4 <- fread("/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.merge.chr.bed") %>% as.data.frame()
mergeG4$length = mergeG4$V3-mergeG4$V2 + 1 ##在处理区间时，通常使用的是左闭右开的区间表示方式
n = c(nrow(mergeG4),mean(mergeG4$length))
mergeG4[mergeG4$length>50,"length"] = 50
ggplot(mergeG4,aes(length)) +
  geom_histogram(binwidth = 1,aes(y=..density..*100)) + #添加直方图，每个直方条的宽度为1个单位，并且将y轴的值转换为密度的百分比形式
  xlab("G4 size (nt)") +
  ylab("Density (%)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 12),axis.text.x = element_text(size = 10,colour = "black"),
    axis.title.y = element_text(size = 12),axis.text.y = element_text(size = 10,colour = "black")) +
  scale_y_continuous(expand = c(0,0),limits = c(0,6))

ggplot(mergeG4,aes(length)) +
  geom_bar(aes(y = (..count..)/sum(..count..)*100),width = 0.75,fill="#F0988C") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 12),axis.text.x = element_text(size = 10,color = "black"),
    axis.title.y = element_text(size = 12),axis.text.y = element_text(size = 10,colour = "black")
  ) +
  scale_y_continuous(expand = c(0,0),limits = c(0,5.5)) +
  xlab("G4 size (nt)") +
  ylab("Density (%)") +
  annotate(geom = "text", x=40, y=4.5,size=4,
           label=paste0("Total eG4","\n",n[1],"\n","Mean size","\n",ceiling(n[2])," nts"))
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"LengthDistributionMergeG4.Dmel.pdf"),
       device = "pdf",width = 3.2,height = 3.2)
#* non_eG4------------------------------------------------------------------------------------
non_eG4 <- fread("/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.non_eG4.chr.bed") %>% as.data.frame()
non_eG4$length = non_eG4$V3-non_eG4$V2 + 1 ##在处理区间时，通常使用的是左闭右开的区间表示方式
n = c(nrow(non_eG4),mean(non_eG4$length))
non_eG4[non_eG4$length>50,"length"] = 50
ggplot(non_eG4,aes(length)) +
  geom_bar(aes(y = (..count..)/sum(..count..)*100),width = 0.75,fill="#A1A9D0") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 12),axis.text.x = element_text(size = 10,color = "black"),
    axis.title.y = element_text(size = 12),axis.text.y = element_text(size = 10,colour = "black")
  ) +
  scale_y_continuous(expand = c(0,0),limits = c(0,5.5)) +
  xlab("G4 size (nt)") +
  ylab("Density (%)") +
  annotate(geom = "text", x=40, y=4.5,size=4,
           label=paste0("Total non-eG4","\n",n[1],"\n","Mean size","\n",ceiling(n[2])," nts"))
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"LengthDistributionnon_eG4.Dmel.pdf"),
       device = "pdf",width = 3.2,height = 3.2)
#* Length Structural stability -----------------------------------------------------
cor <- cor.test(mergeG4$length,mergeG4$V5) ##计算这两个向量之间的皮尔逊相关系数(线性关系)
cor
# 皮尔逊相关系数衡量两个数值变量之间的线性关系强度和方向。它的取值范围在 -1 到 1 之间，其中：
# 1 表示完全正相关，即当一个变量增加时，另一个变量也增加。
# -1 表示完全负相关，即当一个变量增加时，另一个变量减少。
# 0 表示没有线性关系，即两个变量之间没有显著的线性趋势。
cor$estimate
pvalue <- sprintf("%.2e",cor$p.value)
# "%.2f" 和 "%.2e" 在格式化时分别使用浮点数和科学计数法的形式，并控制保留的小数位数或有效数字位数。

ggplot(mergeG4,aes(length,V5))+
  geom_point()+
  geom_smooth(method = "lm") + # 添加线性拟合曲线
  theme_bw() +
  ylab("Structural stability") +
  theme(
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 12),axis.text.x = element_text(size = 10,colour = "black"),
    axis.title.y = element_text(size = 12),axis.text.y = element_text(size = 10,colour = "black")
  ) +
  annotate(geom = "text",x=42,y=166,size=4,
           label=paste0("R=",round(cor$estimate,2),","," P=",pvalue))
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"LengthStructureMergeG4.Dmel.pdf"),
       device = "pdf",width = 5.2,height = 4.2)

#### 长度排序 ####
mergeG4 %<>% arrange(mergeG4$length)
mergeG4 %<>% arrange(length) ##两种都可以
mergeG4$group = cut(mergeG4$length, breaks = c(10,25,28,34,42,135), labels = c("SS","S","M","L","LL"),
                  right = FALSE)
my_comparisons = list(c("SS","S"),c("S","M"),c("M","L"),c("L","LL"))
ggplot(mergeG4,aes(group,V5,fill=group)) +
  geom_boxplot(outlier.colour = "white",notch = TRUE) + ##缺口
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 12),axis.text.x = element_text(size = 10,colour = "black"),
    axis.title.y = element_text(size = 12),axis.text.y = element_text(size = 10,colour = "black")
  ) +
  coord_cartesian(ylim = c(50,120)) +
  ylab("Structural stability") +
  scale_fill_manual(values = c("#c6dbef","#9ecae1","#6baed6","#3182bd","#08519c")) +
  scale_color_manual(values = c("#c6dbef","#9ecae1","#6baed6","#3182bd","#08519c"))+
  stat_compare_means(comparisons = my_comparisons, tip.length = 0,
                     bracket.size = 0.0,label.y = rep(110,4),label = "p.signif",size=4) +
  guides(fill="none")
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"LengthgroupStructureMergeG4.Dmel.pdf"),
       device = "pdf",width = 5.2,height = 4.2)
apply(mergeG4,mergeG4$V5,median) ##箱线图中的中间线即为中位数，它代表数据的中间值
apply(mergeG4,mergeG4$V5,mean) ##看数据怎么变化的，谁高谁低用平均值
##记住！显著性通常是根据数据分布得出的，所以长度和结构稳定性是正相关的

#* Number -----------------------------------------------------
mergeG4$V1 %<>% gsub("chr","",.)
ggplot(mergeG4,aes(V1,fill=V1)) +
  geom_bar(aes(y = ((..count..)/10^4)),position = "dodge") +
  xlab("Chromosome arm") +
  ylab(bquote(Number~(x10^4))) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12)
  ) +
  scale_y_continuous(expand = c(0,0),limits = c(0,0.5)) +
  scale_fill_manual(values = c(rep("#06577A",5),"#d6604d","#06577A")) +
  guides(fill="none")
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"NumberMergeG4.Dmel.pdf"),
       device = "pdf",width = 4.2,height = 3)

#* Density -----------------------------------------------------------------
chr_len <- fread("/home/yuss/flyG4/result/PQS/001.2.chr.length.txt") %>% as.data.frame()
G4number = mergeG4 %>% group_by(V1) %>% summarise(num=n()) %>% as.data.frame() ##summarise(num=n()) 的作用就是在每个数据组内计算行数，将结果存储在一个名为 num 的列中，用于表示每个数据组的行数。
G4number = merge(G4number,chr_len,by.x="V1",by.y="V1")
G4number$density = G4number$num/G4number$V2*10^3
ggplot(G4number,aes(V1,density,fill=V1))+geom_bar(stat = "identity",position = "dodge")+
  xlab("Chromosome arm") + ylab("Density (No./kb)")+
  theme_classic()+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=10,colour = "black"),axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=10,colour = "black"),axis.title.y = element_text(size=12))+
  scale_fill_manual(values=c(rep("#06577A",5),"#d6604d","#06577A"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,0.16))+
  guides(fill="none")
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"DensityMergeG4.Dmel.pdf"),
       device = "pdf",width = 4.2,height = 3)

chr_gc <- fread("/home/yuss/flyG4/result/PQS/001.2.chr.GC.txt") %>% as.data.frame()
G4number$gc <- chr_gc[match(G4number$V1,chr_gc$`#1_usercol`),5]
G4number$norden = G4number$density / G4number$gc
ggplot(G4number,aes(V1,norden,fill=V1))+geom_bar(stat = "identity",position = "dodge")+
  xlab("Chromosome arm") + ylab("Normalized density (No./kb/GC%)")+
  theme_classic()+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=10,colour = "black"),axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=10,colour = "black"),axis.title.y = element_text(size=12))+
  scale_fill_manual(values=c(rep("#06577A",5),"#d6604d","#06577A"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,0.4))+
  guides(fill="none")
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"NorDensityMergeG4.Dmel.pdf"),
       device = "pdf",width = 4.2,height = 3)

#* Structure ----------------------------------------------------------------------------
p = c()
for (i in setdiff(unique(mergeG4$V1),"X")) {
  p[i] = wilcox.test(mergeG4[mergeG4$V1==i,5],mergeG4[mergeG4$V1=="X",5],alternative = "less")$p.value
}
if (all(p<0.001)) {
  ggplot(mergeG4,aes(V1,V5,color=V1)) + geom_boxplot(outlier.colour = "white",notch = TRUE)
} ##做不了，因为p值大，不具有显著性
# p
# 2L           2R           3L           3R            4            Y 
# 6.383692e-02 3.653016e-05 1.721446e-03 3.680628e-04 4.413895e-02 6.332184e-01 
my_comparisons = list(c("2L","2R"),c("2R","3L"),c("3L","3R"),c("3R","4"),c("4","X"),c("X","Y"))
ggplot(mergeG4,aes(V1,V5,color=V1)) + geom_boxplot(outlier.colour = "white",notch = TRUE) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 12),axis.text.x = element_text(size = 10,colour = "black"),
    axis.title.y = element_text(size = 12),axis.text.y = element_text(size = 10,colour = "black")
  ) +
  xlab("Chromosome arm") + ylab("Structural stability") +
  coord_cartesian(ylim = c(50,120)) +
  scale_color_manual(values=c(rep("#06577A",5),"#d6604d","#06577A")) +
  stat_compare_means(comparisons = my_comparisons, tip.length = 0,
                     label.y = rep(c(100,110),3),size=4) +
  guides(color="none")
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,".StructureMergeG4.Dmel.pdf"),
       device = "pdf",width = 4.2,height = 3)

#* Length --------------------------------------------------------
ggplot(mergeG4,aes(V1,length,color=V1)) + geom_boxplot(outlier.colour = "white",notch = TRUE) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 12),axis.text.x = element_text(size = 10,colour = "black"),
    axis.title.y = element_text(size = 12),axis.text.y = element_text(size = 10,colour = "black")
  ) +
  xlab("Chromosome arm") + ylab("Length") +
  coord_cartesian(ylim = c(20,70)) +
  scale_color_manual(values=c(rep("#06577A",5),"#d6604d","#06577A")) +
  stat_compare_means(comparisons = my_comparisons, tip.length = 0,
                     label.y = rep(c(60,65),3),size=4) +
  guides(color="none")
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,".LengthMergeG4.Dmel.pdf"),
       device = "pdf",width = 4.2,height = 3)

由于Y和4号染色体数量太少会影响整体的结果，所以去掉，重新做数量等图
#* Structure after strike out Y、4 Chromosome --------------------------------------------------------
sg4 <- mergeG4[c(mergeG4$V1!="Y" & mergeG4$V1!="4"),]
p = c()
for (i in setdiff(unique(sg4$V1),"X")) {
  p[i] = wilcox.test(sg4[sg4$V1==i,5],sg4[sg4$V1=="X",5],alternative = "less")$p.value
}
if (all(p<0.05)) {
  ggplot(mergeG4,aes(V1,V5,color=V1)) + geom_boxplot(outlier.colour = "white",notch = TRUE)
} ##做不了，因为p值大，不具有显著性
# p
# 2L           2R           3L           3R            4            Y 
# 6.383692e-02 3.653016e-05 1.721446e-03 3.680628e-04 4.413895e-02 6.332184e-01 
my_comparisons = list(c("2L","2R"),c("2R","3L"),c("3L","3R"),c("3R","X"))
ggplot(sg4,aes(V1,V5,color=V1)) + geom_boxplot(outlier.colour = "white",notch = TRUE) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 12),axis.text.x = element_text(size = 10,colour = "black"),
    axis.title.y = element_text(size = 12),axis.text.y = element_text(size = 10,colour = "black")
  ) +
  xlab("Chromosome arm") + ylab("Structural stability") +
  coord_cartesian(ylim = c(50,120)) +
  scale_color_manual(values=c(rep("#06577A",4),"#d6604d")) +
  stat_compare_means(comparisons = my_comparisons, tip.length = 0,
                     label.y = rep(c(100,110),3),size=4) +
  guides(color="none")
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,".StructureMergeG4out.Dmel.pdf"),
       device = "pdf",width = 4.2,height = 3)

#* Length after strike out Y、4 Chromosome --------------------------------------------------------
ggplot(sg4,aes(V1,length,color=V1)) + geom_boxplot(outlier.colour = "white",notch = TRUE) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 12),axis.text.x = element_text(size = 10,colour = "black"),
    axis.title.y = element_text(size = 12),axis.text.y = element_text(size = 10,colour = "black")
  ) +
  xlab("Chromosome arm") + ylab("Length") +
  coord_cartesian(ylim = c(20,60)) +
  scale_color_manual(values=c(rep("#06577A",4),"#d6604d")) +
  stat_compare_means(comparisons = my_comparisons, tip.length = 0,
                     label.y = rep(c(53,57),3),size=4) +
  guides(color="none")
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,".LengthMergeG4out.Dmel.pdf"),
       device = "pdf",width = 4.2,height = 3)
