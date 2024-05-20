#G4-CUT&Tag峰所在基因不同位置的数量占比(在基因组特征区域的分布)
#注释参考文件，即需要一个包含注释信息的TxDb对象
#注释信息一般要包含基因的起始终止，基因的外显子、内含子及它们的起始终止、非编码区域位置、功能元件的位置等
#BiocManager::install("TxDb.Dmelanogaster.UCSC.dm6.ensGene") #同一物种的不同版本TxDb，用哪个取决于上游分析比对使用的哪个版本的基因组
#TxDb注释文件

rm(list = ls());gc();rm(list = ls())
Num = "001.3."
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
Dmel.txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
#peaks注释,tss(转录起始位点)选取上下游2.5Kb，可自行调整，这里默认promoter是上游2.5Kb
anno <- annotatePeak("/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.kc_specific.chr.bed",
                     tssRegion = c(-2000,2000),TxDb = Dmel.txdb)
anno2 <- annotatePeak("/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.s2_specific.chr.bed",
                      tssRegion = c(-2000,2000),TxDb = Dmel.txdb)
anno3 <- annotatePeak("/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.merge.chr.bed",
                      tssRegion = c(-2000,2000),TxDb = Dmel.txdb)
anno4 <- annotatePeak("/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.overlap.chr.bed",
                      tssRegion = c(-2000,2000),TxDb = Dmel.txdb)
anno5 <- annotatePeak("/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.non_eG4.chr.bed",
                      tssRegion = c(-2000,2000),TxDb = Dmel.txdb)
df1 <- anno@annoStat %>% as.data.frame()
df2 <- anno2@annoStat %>% as.data.frame()
df3 <- anno3@annoStat %>% as.data.frame()
df4 <- anno4@annoStat %>% as.data.frame()
df5 <- anno5@annoStat %>% as.data.frame()

df1$Feature %<>% as.character()
df1[nrow(df1)+1,] = c("Exon",df1[df1$Feature=="1st Exon",2]+df1[df1$Feature=="Other Exon",2])
df1[nrow(df1)+1,] = c("Intron", as.numeric(df1[df1$Feature=="1st Intron",2])+as.numeric(df1[df1$Feature=="Other Intron",2]))
df1[nrow(df1)+1,] = c("Promoter",as.numeric(df1[df1$Feature=="Promoter (<=1kb)",2])+as.numeric(df1[df1$Feature=="Promoter (1-2kb)",2]))
df1 %<>% dplyr::filter(Feature %in% c("Promoter","5' UTR","3' UTR","Exon","Intron","Distal Intergenic","Downstream (<=300)"))
df1$Frequency %<>% as.numeric()
df1
df1$Feature %<>% gsub("Downstream (<=300)","Downstream",.,fixed = TRUE) 
df1$Feature <- factor(df1$Feature,levels = c("Downstream","5' UTR","3' UTR","Exon","Distal Intergenic","Intron","Promoter"))
df1

df2$Feature %<>% as.character()
#df2[nrow(df2)+1,] = c("Exon",df2[df2$Feature=="1st Exon",2]+df2[df2$Feature=="Other Exon",2])
df2[nrow(df2)+1,] = c("Intron", as.numeric(df2[df2$Feature=="1st Intron",2])+as.numeric(df2[df2$Feature=="Other Intron",2]))
df2[nrow(df2)+1,] = c("Promoter",as.numeric(df2[df2$Feature=="Promoter (<=1kb)",2])+as.numeric(df2[df2$Feature=="Promoter (1-2kb)",2]))
df2$Feature %<>% gsub("Other Exon","Exon",.,fixed = TRUE) 
df2 %<>% dplyr::filter(Feature %in% c("Promoter","5' UTR","3' UTR","Exon","Intron","Distal Intergenic","Downstream (<=300)"))
df2$Frequency %<>% as.numeric()
df2
df2$Feature %<>% gsub("Downstream (<=300)","Downstream",.,fixed = TRUE) 
df2$Feature <- factor(df2$Feature,levels = c("Downstream","5' UTR","3' UTR","Exon","Distal Intergenic","Intron","Promoter"))
df2

df3$Feature %<>% as.character()
df3[nrow(df3)+1,] = c("Exon",df3[df3$Feature=="1st Exon",2]+df3[df3$Feature=="Other Exon",2])
df3[nrow(df3)+1,] = c("Intron", as.numeric(df3[df3$Feature=="1st Intron",2])+as.numeric(df3[df3$Feature=="Other Intron",2]))
df3[nrow(df3)+1,] = c("Promoter",as.numeric(df3[df3$Feature=="Promoter (<=1kb)",2])+as.numeric(df3[df3$Feature=="Promoter (1-2kb)",2]))
df3 %<>% dplyr::filter(Feature %in% c("Promoter","5' UTR","3' UTR","Exon","Intron","Distal Intergenic","Downstream (<=300)"))
df3$Frequency %<>% as.numeric()
df3
df3$Feature %<>% gsub("Downstream (<=300)","Downstream",.,fixed = TRUE) 
df3$Feature <- factor(df3$Feature,levels = c("Downstream","5' UTR","3' UTR","Exon","Distal Intergenic","Intron","Promoter"))
df3

df4$Feature %<>% as.character()
df4[nrow(df4)+1,] = c("Exon",df4[df4$Feature=="1st Exon",2]+df4[df4$Feature=="Other Exon",2])
df4[nrow(df4)+1,] = c("Intron", as.numeric(df4[df4$Feature=="1st Intron",2])+as.numeric(df4[df4$Feature=="Other Intron",2]))
df4[nrow(df4)+1,] = c("Promoter",as.numeric(df4[df4$Feature=="Promoter (<=1kb)",2])+as.numeric(df4[df4$Feature=="Promoter (1-2kb)",2]))
df4 %<>% dplyr::filter(Feature %in% c("Promoter","5' UTR","3' UTR","Exon","Intron","Distal Intergenic","Downstream (<=300)"))
df4$Frequency %<>% as.numeric()
df4
df4$Feature %<>% gsub("Downstream (<=300)","Downstream",.,fixed = TRUE) 
df4$Feature <- factor(df4$Feature,levels = c("Downstream","5' UTR","3' UTR","Exon","Distal Intergenic","Intron","Promoter"))
df4

df5$Feature %<>% as.character()
df5[nrow(df5)+1,] = c("Exon",df5[df5$Feature=="1st Exon",2]+df5[df5$Feature=="Other Exon",2])
df5[nrow(df5)+1,] = c("Intron", as.numeric(df5[df5$Feature=="1st Intron",2])+as.numeric(df5[df5$Feature=="Other Intron",2]))
df5[nrow(df5)+1,] = c("Promoter",as.numeric(df5[df5$Feature=="Promoter (<=1kb)",2])+as.numeric(df5[df5$Feature=="Promoter (1-2kb)",2]))
df5 %<>% dplyr::filter(Feature %in% c("Promoter","5' UTR","3' UTR","Exon","Intron","Distal Intergenic","Downstream (<=300)"))
df5$Frequency %<>% as.numeric()
df5
df5$Feature %<>% gsub("Downstream (<=300)","Downstream",.,fixed = TRUE) 
df5$Feature <- factor(df5$Feature,levels = c("Downstream","5' UTR","3' UTR","Exon","Distal Intergenic","Intron","Promoter"))
df5

df1$Type <- 'kc_specific'
df2$Type <- 's2_specific'
df3$Type <- 'merge'
df4$Type <- 'overlap'
df5$Type <- 'non_eG4'

df <- rbind(df1,df2,df3,df4,df5)
df$Type %<>% factor(.,levels = c("non_eG4","kc_specific","s2_specific","overlap","merge"))

b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(5)
color <- c('#4d4d4d',b)
#画图
library(ggplot2)
ggplot(df,aes(x = Feature,y= Frequency,fill = Type)) +
  geom_bar(position = 'dodge',stat = 'identity',color = 'white') +
  scale_fill_manual(values = color)+
  scale_y_continuous(expand = c(0,0)) +
  cowplot::theme_half_open() +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_blank(),
        axis.text.x = element_text(vjust = 0.5,hjust = 1,size = 14),
        axis.text.y = element_text(size=14),
        axis.ticks.x = element_blank(), #修改x轴刻度线
        legend.position = 'top',
        legend.title=element_blank(),
        plot.margin = margin(t = 10,  # 顶部边缘距离
                             r = 10,  # 右边边缘距离
                             b = 10,  # 底部边缘距离
                             l = 10)) +
  ylab("Frequency (%)") +
  geom_text(aes(label=str_c(round(Frequency,2),"%"),y=Frequency + 2), color="black", size=3.5,
            position = position_dodge(.9)) + #str_c函数是将多个字符串合并成一个字符串，y=Frequency + 3是标签高低位置设置，position = position_dodge(.9)是标签左右位置设置
  coord_flip(ylim = c(0,50)) #x 和 y 翻转的笛卡尔坐标
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"PQS.GenomicAnnotation.pdf"),
        device = "pdf",width = 10,height = 6.5)

path <- "/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed"
files <- dir(path)
filepath <- sapply(files, function(x)
  {paste(path,x,sep = '/')
})
peak=list(non_eG4=filepath[[3]],kc_specific=filepath[[1]],s2_specific=filepath[[5]],overlap=filepath[[4]],merge=filepath[[2]])
promoter <- getPromoters(TxDb = Dmel.txdb,upstream = 2000,downstream = 2000)
tagMatrixList <- lapply(peak, getTagMatrix, windows=promoter)
b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c('#4d4d4d',b)

p <- plotAvgProf(tagMatrixList, xlim=c(-2000, 2000),
                 xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency") +
  theme_bw()+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size= 16),
        axis.text = element_text(size= 14),
        legend.position=c(0.1,0.8),legend.title = element_blank(),
        legend.background = element_rect(fill = NA)) +labs(color = 'Group')+
  scale_fill_manual(values = color,aesthetics = "color") ##折线图

rm(list = ls());gc();rm(list = ls())
Num = "001.3."
##awk '{print "chr"$0}' "/home/yuss/flyG4/result/PQS/5.type.bed/001.2.s2_all.bed" > /home/yuss/flyG4/result/PQS/5.type.bed/001.2.s2_all.chr.bed
#### genome annotation (s2 all eG4和kc all eG4)####
#*s2 all eG4基因组注释----------------------------------------------------------------------------------
##创建转录组数据库（TxDb)
##直接使用TxDb包中提供的函数来构建TxDb
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
Dmel.txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
##也可以用makeTxDbFromGFF()函数从GFF文件中创建TxDb   Dmel.txdb <-makeTxDbFromGFF("/home/qians/Quadruplex/Input/Ref/Fly/dmel-all-r6.19.gff",format = "gff3")
library(ChIPseeker)
s2all.anno <- annotatePeak("/home/yuss/flyG4/result/PQS/5.type.bed/001.2.s2_all.chr.bed",
                           tssRegion = c(-2000,2000),TxDb = Dmel.txdb)
s2all.anno
## 可视化
plotAnnoBar(s2all.anno)
plotAnnoPie(s2all.anno)
vennpie(s2all.anno)

df = s2all.anno@annoStat
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

df$lab <- paste0(round(df$Frequency,2),"%")
ggdonutchart(df, "Frequency",
             label = "lab",
             lab.pos = "in",
             fill = "Feature",
             color = "white")


library(ggplot2)
library(ggforce)
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.direction = "horizontal",
        legend.position = c(0.5, 0)
  )+ #去除没用的ggplot背景
  guides(fill = guide_legend(reverse = TRUE))+ #使用guides语法反转图例的项目顺序
  xlab("")+ylab('')+#坐标轴
  scale_fill_manual(values = c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', 
                               '#D6E7A3', '#57C3F3','#476D87'))+ #添加颜色
  geom_arc_bar(data=df,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=Frequency,fill=Feature),
                   color = "transparent"
  )+#饼图（r0=0）
  annotate("text",x=1.95,y=1.0,label="38.09%",angle=-62)+
  annotate("text",x=0,y=-2.2,label="28.82%",angle=0)+
  annotate("text",x=-2.2,y=0.3,label="21.51%",angle=80)+
  annotate("text",x=-0.9,y=2,label="9.31%",angle=20)+
  annotate("text",x=-0.1,y=2.18,label="1.78%",angle=5)
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"LocationS2allG4.Dmel.pdf"),
       device = "pdf",width = 9,height = 9.2)

#*kc all eG4基因组注释----------------------------------------------------------------------------------
kcall.anno <- annotatePeak("/home/yuss/flyG4/result/PQS/5.type.bed/001.2.kc_all.chr.bed",
                           tssRegion = c(-2000,2000),TxDb = Dmel.txdb)
kcall.anno
## 可视化
plotAnnoBar(kcall.anno)
plotAnnoPie(kcall.anno)
vennpie(kcall.anno)

df = kcall.anno@annoStat
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
df$lab <- paste0(round(df$Frequency,2),"%")
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.direction = "horizontal",
        legend.position = c(0.5, 0)
        )+ #去除没用的ggplot背景
  guides(fill = guide_legend(reverse = TRUE))+ #使用guides语法反转图例的项目顺序
  xlab("")+ylab('')+#坐标轴
  scale_fill_manual(values = c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', 
                               '#D6E7A3', '#57C3F3','#476D87'))+ #添加颜色
  geom_arc_bar(data=df,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=Frequency,fill=Feature),
               color = "transparent"
  )+#饼图（r0=0）
  annotate("text",x=1.95,y=1.0,label="36.06%",angle=-62)+
  annotate("text",x=0,y=-2.2,label="30.18%",angle=0)+
  annotate("text",x=-2.2,y=0.3,label="22.69%",angle=80)+
  annotate("text",x=-0.9,y=2,label="8.76%",angle=20)+
  annotate("text",x=-0.1,y=2.18,label="1.75%",angle=5)
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"LocationKcallG4.Dmel.pdf"),
device = "pdf",width = 9,height = 9.2)


#*kc specific eG4基因组注释----------------------------------------------------------------------------------
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
Dmel.txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
kcspe.anno <- annotatePeak("/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.kc_specific.chr.bed",
                           tssRegion = c(-2000,2000),TxDb = Dmel.txdb)
kcspe.anno
## 可视化
plotAnnoBar(kcspe.anno)
plotAnnoPie(kcspe.anno)
vennpie(kcspe.anno)

df = kcspe.anno@annoStat
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
df$lab <- paste0(round(df$Frequency,2),"%")
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.direction = "horizontal",
        legend.position = c(0.5, 0)
  )+ #去除没用的ggplot背景
  guides(fill = guide_legend(reverse = TRUE))+ #使用guides语法反转图例的项目顺序
  xlab("")+ylab('')+#坐标轴
  scale_fill_manual(values = c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', 
                               '#D6E7A3', '#57C3F3','#476D87'))+ #添加颜色
  geom_arc_bar(data=df,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=Frequency,fill=Feature),
               color = "transparent"
  )+#饼图（r0=0）
  annotate("text",x=1.95,y=1.0,label="36.06%",angle=-62)+
  annotate("text",x=0,y=-2.2,label="30.18%",angle=0)+
  annotate("text",x=-2.2,y=0.3,label="22.69%",angle=80)+
  annotate("text",x=-0.9,y=2,label="8.76%",angle=20)+
  annotate("text",x=-0.1,y=2.18,label="1.75%",angle=5)
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"LocationKcallG4.Dmel.pdf"),
       device = "pdf",width = 9,height = 9.2)

s2spe.anno <- annotatePeak("/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.s2_specific.chr.bed",
                           tssRegion = c(-2000,2000),TxDb = Dmel.txdb)

s2spe.anno
## 可视化
plotAnnoBar(s2all.anno)
plotAnnoPie(s2all.anno)
vennpie(s2all.anno)

df = s2spe.anno@annoStat
class(df$Feature) ##factor
df$Feature %<>% as.character()
df[nrow(df)+1,] = c("Exon",df[df$Feature=="Other Exon",2])
df[nrow(df)+1,] = c("Intron", as.numeric(df[df$Feature=="1st Intron",2])+as.numeric(df[df$Feature=="Other Intron",2]))
df[nrow(df)+1,] = c("Promoter",as.numeric(df[df$Feature=="Promoter (<=1kb)",2])+as.numeric(df[df$Feature=="Promoter (1-2kb)",2]))
df %<>% dplyr::filter(Feature %in% c("Promoter","5' UTR","3' UTR","Exon","Intron","Distal Intergenic","Downstream (<=300)"))
df$Frequency %<>% as.numeric()
df
df %<>% arrange(-Frequency)
## df 中的 Feature 列中的特定字符串 "Downstream (<=300)" 替换为 "Downstream"，并将修改后的结果存回 df$Feature 列,fixed = TRUE 参数表示使用字符串匹配而不是正则表达式匹配,接着把这列设置为因子，rev(.) 表示将列的因子水平逆序排列。
df$Feature %<>% gsub("Downstream (<=300)","Downstream",.,fixed = TRUE) %>% factor(.,levels = rev(.))
library(ggplot2)
library(ggforce)
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.direction = "horizontal",
        legend.position = c(0.5, 0)
  )+ #去除没用的ggplot背景
  guides(fill = guide_legend(reverse = TRUE))+ #使用guides语法反转图例的项目顺序
  xlab("")+ylab('')+#坐标轴
  scale_fill_manual(values = c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', 
                               '#D6E7A3', '#57C3F3','#476D87'))+ #添加颜色
  geom_arc_bar(data=df,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=Frequency,fill=Feature),
               color = "transparent"
  )+#饼图（r0=0）
  annotate("text",x=1.95,y=1.0,label="38.09%",angle=-62)+
  annotate("text",x=0,y=-2.2,label="28.82%",angle=0)+
  annotate("text",x=-2.2,y=0.3,label="21.51%",angle=80)+
  annotate("text",x=-0.9,y=2,label="9.31%",angle=20)+
  annotate("text",x=-0.1,y=2.18,label="1.78%",angle=5)
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"LocationS2allG4.Dmel.pdf"),
       device = "pdf",width = 9,height = 9.2)
