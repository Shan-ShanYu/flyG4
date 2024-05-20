rm(list = ls());gc();rm(list = ls())#清空
Num = "001.2."
setwd("/home/yuss/flyG4/result/PQS")
#### pqs和peak交集后G4分类####
#读取pqs和peak交集的表
pqs_kc_1 <- fread('/home/yuss/flyG4/result/PQS/001.2.pqs.KCrep1_R1_peaks.bed') %>% data.frame()
pqs_kc_2 <- fread('/home/yuss/flyG4/result/PQS/001.2.pqs.KCrep2_R1_peaks.bed') %>% data.frame()
pqs_s2_1 <- fread('/home/yuss/flyG4/result/PQS/001.2.pqs.S2rep1_R1_peaks.bed') %>% data.frame()
pqs_s2_2 <- fread('/home/yuss/flyG4/result/PQS/001.2.pqs.S2rep2_R1_peaks.bed') %>% data.frame()

#合并交集（pqs信息在前，交集信息在后）
library(dplyr)
pqs_peak <- bind_cols(pqs_kc_1,pqs_kc_2$V7,pqs_s2_1$V7,pqs_s2_2$V7)

#修改列名
colnames(pqs_peak) <- c("chr","start","end","id","score","strand","kc_1","kc_2","s2_1","s2_2")
##改单个列名colnames(pqs_peak)[7] <- 'kc_1'

#合并两个重复
pqs_peak$kc <- ifelse(rowSums(pqs_peak[,c("kc_1", "kc_2")])==2,1,0)
pqs_peak$s2 <- ifelse(rowSums(pqs_peak[,c("s2_1", "s2_2")])==2,1,0)

#改$s2
pqs_peak$s2 <- ifelse(pqs_peak$s2==1,2,0)
##修改为数字类型pqs_peak$s2 <- as.numeric(pqs_peak$s2)

#求kc和s2的sum
pqs_peak$sum <- rowSums(pqs_peak[,11:12])
kc_specific <- pqs_peak[pqs_peak$sum==1,]
s2_specific <- pqs_peak[pqs_peak$sum==2,]
overlap <- pqs_peak[pqs_peak$sum==3,]
merge <- pqs_peak[pqs_peak$sum!=0,]
non_eG4 <-pqs_peak[pqs_peak$sum==0,]
kc_all <- pqs_peak[pqs_peak$kc!=0,]
s2_all <- pqs_peak[pqs_peak$s2!=0,]
#保存
write.table(pqs_peak,file = paste0("/home/yuss/flyG4/result/PQS/",Num,"pqs_peak.bed"),sep = '\t',
            col.names = T,row.names = F,quote = F)
write.table(kc_specific,file = paste0("/home/yuss/flyG4/result/PQS/",Num,"kc_specific.bed"),sep = '\t',
            col.names = T,row.names = F,quote = F)
write.table(s2_specific,file = paste0("/home/yuss/flyG4/result/PQS/",Num,"s2_specific.bed"),sep = '\t',
            col.names = T,row.names = F,quote = F)
write.table(overlap,file = paste0("/home/yuss/flyG4/result/PQS/",Num,"overlap.bed"),sep = '\t',
            col.names = T,row.names = F,quote = F)
write.table(merge,file = paste0("/home/yuss/flyG4/result/PQS/",Num,"merge.bed"),sep = '\t',
            col.names = T,row.names = F,quote = F)
write.table(non_eG4,file = paste0("/home/yuss/flyG4/result/PQS/",Num,"non_eG4.bed"),sep = '\t',
            col.names = T,row.names = F,quote = F)
write.table(kc_all,file = paste0("/home/yuss/flyG4/result/PQS/",Num,"kc_all.bed"),sep = '\t',
            col.names = T,row.names = F,quote = F)
write.table(s2_all,file = paste0("/home/yuss/flyG4/result/PQS/",Num,"s2_all.bed"),sep = '\t',
            col.names = T,row.names = F,quote = F)

#按行合并non_eG4,kc_specific,s2_specific,overlap,merge
kc_specific$type <- "kc_specific"
s2_specific$type <- "s2_specific"
overlap$type <- "overlap" 
merge$type <- "merge" 
non_eG4$type <- "non_eG4"
merge.all <- bind_rows(non_eG4,kc_specific,s2_specific,overlap,merge) ##bind_rows()按行合并，把表竖起来
write.table(merge.all,file = paste0("/home/yuss/flyG4/result/PQS/",Num,"merge.all.bed"),sep = '\t',
            col.names = T,row.names = F,quote = F)
#### G4的数量 ####
#*每类G4的数量---------------------------------------------------------------------------------
merge.all <- fread("/home/yuss/flyG4/result/PQS/001.2.merge.all.bed") %>% as.data.frame()
#计数
num <- table(merge.all['type']) %>% as.data.frame() ##table()来计算mou一列中的出现次数
num$type <- factor(num$type,levels =c("non_eG4","kc_specific","s2_specific","overlap","merge"))

#non_eG4,kc_specific,s2_specific,overlap,merge的数量图
library(ggplot2)
library(ggsci)
b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c('#4d4d4d',b)
ggplot(num, aes(x=type,y=Freq,fill=type)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(0.7),width = 0.7,color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值
  coord_cartesian(ylim = c(0,30000)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = color) +
  geom_text(aes(label=Freq,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题(左下边框，没有网格线)
  ylab("Number of G4") + ##通过bquote函数给图标签添加上下标
  theme(axis.title.y = element_text(size = 16), ##y坐标轴标题字体大小
        axis.title.x = element_blank(), ##删除x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.position = "none") 
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"G4Number.pdf"),
       device = "pdf",width = 6.2, height = 4)    
# ggplot(data = b,aes(x = Var1,y= num,fill = Var1)) +
#   geom_bar(stat = 'identity',position=position_dodge(0.7),width = 0.7,color = 'white') +
#   coord_cartesian(ylim = c(0,1.3)) + #坐标轴范围
#   scale_y_continuous(expand = c(0,0)) +
#   scale_fill_manual(values = a) + #用于面的填充色
#   geom_text(aes(label=num,vjust = -0.5), color="black", size=4) + 
#   cowplot::theme_half_open()+
#   ylab(bquote(Number~of~G4~(x10^5))) + 
#   theme(axis.title.y = element_text(size=16),
#         axis.title.x = element_blank(),
#         axis.text = element_text(size=14),
#         legend.position = "none")

#*non_eG4和mergeG4的数量-----------------------------------------------------------------
num.noneG4.mergeG4 <- num[c(3,2),]
ggplot(num.noneG4.mergeG4, aes(x=type,y=Freq,fill=type)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(0.7),width = 0.7,color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值
  coord_cartesian(ylim = c(0,30000)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("#A1A9D0","#F0988C")) +
  geom_text(aes(label=Freq,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题(左下边框，没有网格线)
  ylab("Number of G4") + ##通过bquote函数给图标签添加上下标
  theme(axis.title.y = element_text(size = 16), ##y坐标轴标题字体大小
        axis.title.x = element_blank(), ##删除x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.position = "none") +
  coord_cartesian(ylim = c(0,30000)) +
  scale_x_discrete(labels = c("non-eG4","eG4"))
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"non_eG4.eG4.Number.pdf"),
       device = "pdf",width = 3.2, height = 3.5)  

#*overlap和细胞特异性G4的数量-----------------------------------------------------------------
num.overlap <- num[c(2:4),]
num.overlap$type <- factor(num.overlap$type,levels = c("overlap","kc_specific","s2_specific"))
#b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)反过来
b = colorRampPalette(colors = c('#C18383', '#FDDBC7'))(3)
ggplot(num.overlap, aes(x=type,y=Freq,fill=type)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(0.7),width = 0.7,color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值
  coord_cartesian(ylim = c(0,12000)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = b) +
  geom_text(aes(label=Freq,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题(左下边框，没有网格线)
  ylab("Number") + ##通过bquote函数给图标签添加上下标
  theme(axis.title.y = element_text(size = 16), ##y坐标轴标题字体大小
        axis.title.x = element_blank(), ##删除x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.position = "none") +
  scale_x_discrete(labels = c("Overlap","Kc specific","S2 specific"))
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"OverlapAndSpecific.Number.pdf"),
       device = "pdf",width = 4.4, height = 4)  

#### G4的长度（end-start） ####
#计算length
merge.all <- fread("/home/yuss/flyG4/result/PQS/001.2.merge.all.bed") %>% as.data.frame()
merge.all$length <- merge.all$end-merge.all$start
merge.all$type<-factor(merge.all$type,levels =c("non_eG4","kc_specific","s2_specific","overlap","merge"))
#画图length
my_comparisons = list(c("non_eG4","kc_specific"),c("kc_specific","s2_specific"),
                      c("s2_specific","overlap"),c("overlap","merge"))
b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c('#4d4d4d',b)
median.length <- aggregate(length ~  type, merge.all, median)
library(ggpubr)
#*每类G4的长度--------------------------------------------------------------------
ggplot(data = merge.all,aes(x=type,y=length,fill=type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,label.y = c(rep(c(60,70),2)),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = color) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Length (bp)") +
  coord_cartesian(ylim = c(0,80)) +
  geom_hline(yintercept = median.length[median.length$type == "non_eG4",2],linetype=2)
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"G4Length.pdf"),
       device = "pdf",width = 6, height = 4)    

#*non_eG4和mergeG4的长度-----------------------------------------------------------------
merge.all.2 <- merge.all[merge.all$type=="non_eG4"|merge.all$type=="merge",] %>% as.data.frame()
my_comparisons = list(c("non_eG4","merge"))
ggplot(data = merge.all.2,aes(x=type,y=length,fill=type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,label.y = 55,
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#A1A9D0","#F0988C")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Length (bp)") +
  coord_cartesian(ylim = c(0,60)) +
  scale_x_discrete(labels = c("non_eG4", "eG4"))
  #geom_hline(yintercept = median.length[median.length$type == "non_eG4",2],linetype=2)
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"non_eG4.eG4Length.pdf"),
       device = "pdf",width = 3.5, height = 3.5) 

#*overlap和细胞特异性G4的长度-----------------------------------------------------------------
merge.all.3 <- merge.all[merge.all$type=="kc_specific"|merge.all$type=="s2_specific"|merge.all$type=="overlap",]
merge.all.3$type <- factor(merge.all.3$type,c("overlap","kc_specific","s2_specific"))
my_comparisons = list(c("overlap","kc_specific"),c("kc_specific","s2_specific"))
#b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)反过来
b = colorRampPalette(colors = c('#C18383', '#FDDBC7'))(3)
median.length <- aggregate(length ~  type, merge.all.3, median)
ggplot(data = merge.all.3,aes(x=type,y=length,fill=type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,label.y = c(51,56),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = b) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Length (bp)") +
  coord_cartesian(ylim = c(5,60)) +
  scale_x_discrete(labels = c("Overlap","Kc specific","S2 specific"))
  #geom_hline(yintercept = median.length[median.length$type == "overlap",2],linetype=2)
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"OverlapAndSpecificG4Length.pdf"),
       device = "pdf",width = 4, height = 3.5) 

#### 每类G4的结构稳定性（score） ####
merge.all <- fread("/home/yuss/flyG4/result/PQS/001.2.merge.all.bed") %>% as.data.frame()
merge.all$type <- factor(merge.all$type,levels =c("non_eG4","kc_specific","s2_specific","overlap","merge"))
median.score <- aggregate(score ~  type, merge.all, median)
library(ggpubr)
ggplot(data = merge.all,aes(x=type,y=score,fill=type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,label.y = c(rep(c(120,135),2)),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4) +
  scale_fill_manual(values = color) +
  cowplot::theme_half_open() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Structural stability") +
  coord_cartesian(ylim = c(50,150)) +
  geom_hline(yintercept = median.score[median.score$type == "non_eG4",2],linetype=2)
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"G4Structure.pdf"),
       device = "pdf",width = 6,height = 4)

tapply(merge.all$score, merge.all$type, median)
table(merge.all$type)
median(merge.all$score)


rm(list = ls());gc();rm(list = ls())#清空
Num = "001.2."
#### 每类G4的在染色体的分布(分组柱形图) ####
merge.all <- fread("/home/yuss/flyG4/result/PQS/001.2.merge.all.bed") %>% as.data.frame()
chr.num <- table(merge.all$chr,merge.all$type)%>% data.frame()
colnames(chr.num) <- c("chr","type","Freq")
chr.num$chr <- factor(chr.num$chr,levels =c('2L','2R','3L','3R',4,'X','Y'))
chr.num$type <- factor(chr.num$type,levels = c("non_eG4","kc_specific","s2_specific","overlap","merge"))

b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c('#4d4d4d',b)
ggplot(chr.num, aes(x=chr,y=Freq,fill=type)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(),width = 0.9,color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值，#width设置矩形条的宽度,color边框颜色
  coord_cartesian(ylim = c(0,8800)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) + ##消除x轴与绘图区的间隙
  scale_fill_manual(values = color) +
  #geom_text(aes(label=Freq,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题
  ylab("Number") +
  xlab("Chromosome arm") +
  #scale_fill_discrete(name="group") + ##修改图例标题名称
  theme(axis.title.y = element_text(size = 16), ##y坐标轴标题字体大小
        axis.title.x = element_text(size = 16), ##x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.position = c(0, 0.97), ##图例在绘图区域的位置 
        legend.direction = "horizontal", ##设置图例水平放置
        legend.title = element_blank()) 
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"G4.ChrNumber.pdf"),
       device = "pdf",width = 6.5, height = 4.2)  
#### 在X染色体上分布的百分比 ####
#*每类G4-----------------------------------------------------------------------------------
##求每类G4在X染色体的数量
x.chr.num <- chr.num[chr.num$chr=="X",]
##分组求和,求每类G4在不同染色体的总数量
chr.num$type <- factor(chr.num$type,c("kc_specific","merge","non_eG4","overlap","s2_specific"))
x.chr.num$sum <- aggregate(chr.num$Freq, by=list(type=chr.num$type),sum)
##求每类G4在X染色体的比例：每类G4在X染色体的数量除以在不同染色体的总数量
x.chr.num$ratio <- x.chr.num$Freq/x.chr.num$sum$x
x.chr.num$percent <- round(x.chr.num$ratio*100, 2)
x.chr.num$type <- factor(x.chr.num$type, levels = c("overlap","non_eG4","kc_specific","s2_specific","merge"))
ggplot(x.chr.num, aes(x=type,y=percent,fill=type)) +
  geom_bar(stat = 'identity',position = position_dodge(),width = 0.7,color = 'white')+
  coord_cartesian(ylim = c(0,50)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("#C18383","#4d4d4d","#FDDBC7","#DFAFA5","#A35762")) +
  geom_text(aes(label=percent,vjust=-0.5),color="black",size=4) +
  cowplot::theme_half_open() +
  ylab("(%)Proportion on the X chromosome")+
  xlab("")+
  theme(axis.title.y = element_text(size = 16),
        axis.text = element_text(size=14),
        legend.position = c(0, 0.97),
        legend.direction = "horizontal",
        legend.title = element_blank())
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"G4.XchrProportion.pdf"),
       device = "pdf",width = 6.4, height = 4.7)

#*overlap和细胞特异性G4-----------------------------------------------------------------------------------
b = colorRampPalette(colors = c('#C18383', '#FDDBC7'))(3)
overlap.x.chr.num <- x.chr.num[c(1,4:5),]
overlap.x.chr.num$type <- factor(overlap.x.chr.num$type,levels = c("overlap","kc_specific","s2_specific"))
ggplot(overlap.x.chr.num, aes(x=type,y=percent,fill=type)) +
  geom_bar(stat = 'identity',position = position_dodge(),width = 0.7,color = 'white')+
  coord_cartesian(ylim = c(0,50)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = b) +
  geom_text(aes(label=percent,vjust=-0.5),color="black",size=4) +
  cowplot::theme_half_open() +
  ylab("(%)Proportion on the X chromosome")+
  xlab("")+
  theme(axis.title.y = element_text(size = 16),
        axis.text = element_text(size=14),
        legend.position = "none")+
  scale_x_discrete(labels = c("Overlap","Kc specific","S2 specific"))
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"OverlapAndSpecificG4.XchrProportion.pdf"),
       device = "pdf",width = 4.4, height = 4.4)

##fisher test 看是否有显著性
x.chr.num$num.chrA <- x.chr.num$sum$x -  x.chr.num$Freq
#num chrX
x.chr.num$Freq
#num chrA
x.chr.num$num.chrA
x.chr.num$type

data1 <- matrix(c(1577,2821,972,8611), nrow = 2)
colnames(data1) <- c("kc","overlap")
rownames(data1) <- c("chr X","chr A")
fisher.test(data1) ##p-value < 2.2e-16

data2 <- matrix(c(231,2190,972,8611), nrow = 2)
colnames(data2) <- c("s2","overlap")
rownames(data2) <- c("chr X","chr A")
fisher.test(data2) ##p-value = 0.4047

data1 <- matrix(c(1577,2821,1203,10801), nrow = 2)
colnames(data1) <- c("kc","other")
rownames(data1) <- c("chr X","chr A")
fisher.test(data1) ##p-value < 2.2e-16

data2 <- matrix(c(231,2190,2549,11432), nrow = 2)
colnames(data2) <- c("s2","other")
rownames(data2) <- c("chr X","chr A")
fisher.test(data2) ##p-value < 2.2e-16

#### Kc.all,S2.all在染色体上的分布 ####
kc.all <- fread("/home/yuss/flyG4/result/PQS/001.2.kc_all.bed") %>% as.data.frame()
s2.all <- fread("/home/yuss/flyG4/result/PQS/001.2.s2_all.bed") %>% as.data.frame()
kc.all$length <- kc.all$end-kc.all$start
s2.all$length <- s2.all$end-s2.all$start
kc.num <- table(kc.all$chr) %>% as.data.frame()
s2.num <- table(s2.all$chr) %>% as.data.frame()
num <- kc.num[c(1:6),]
num$s2 <- s2.num$Freq
colnames(num) <- c("chr","kc","s2")

ggplot(num,aes(x=chr,y=kc))+
  geom_bar(stat = "identity",position = "dodge",fill="#D76364")+
  cowplot::theme_half_open() +
  ylab("Number")+
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=10,colour = "black",angle = 45,hjust = 1,vjust = 1),axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,colour = "black"),axis.title.y = element_text(size=12),
        legend.position = "none")+
  scale_y_continuous(expand = c(0,0),limits = c(0,4000)) +
  labs(title="Kc167 cells(13981)") 
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"KcG4.ChrNumber.pdf"),
       device = "pdf",width = 3.2, height = 3.6)  

ggplot(num,aes(x=chr,y=s2))+
  geom_bar(stat = "identity",position = "dodge",fill="#5F97D3")+
  cowplot::theme_half_open() +
  ylab("Number")+
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=10,colour = "black",angle = 45,hjust = 1,vjust = 1),axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,colour = "black"),axis.title.y = element_text(size=12),
        legend.position = "none")+
  scale_y_continuous(expand = c(0,0),limits = c(0,4000)) +
  labs(title="S2 cells(12004)") #+
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"S2G4.ChrNumber.pdf"),
       device = "pdf",width = 3.2, height = 3.6)

median.length <- aggregate(length ~  chr, kc.all, median)
median.length <- aggregate(length ~  chr, s2.all, median)
kc.all.exclude <- kc.all[kc.all$chr!="Y",]
median.score <- aggregate(score ~  chr, kc.all.exclude, median)
library(ggpubr)
ggplot(data = kc.all.exclude,aes(x=chr,y=score,colour=chr)) +
  geom_boxplot(notch = FALSE,outlier.colour = "white") +
  # stat_compare_means(comparisons = my_comparisons,label.y = c(rep(c(120,135),2)),
  #                    aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4) +
  # scale_fill_manual(values = color) +
  cowplot::theme_half_open() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Structural stability") +
  coord_cartesian(ylim = c(40,110)) +
  labs(title="Kc167 cells")
  geom_hline(yintercept = median.score[median.score$type == "non_eG4",2],linetype=2)
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"KcG4.ChrStability.pdf"),
       device = "pdf",width = 6,height = 4)

median.score <- aggregate(score ~  chr, s2.all, median)
ggplot(data = s2.all,aes(x=chr,y=score,colour=chr)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  # stat_compare_means(comparisons = my_comparisons,label.y = c(rep(c(120,135),2)),
  #                    aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4) +
  # scale_fill_manual(values = color) +
  cowplot::theme_half_open() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Structural stability") +
  coord_cartesian(ylim = c(40,110))# +
geom_hline(yintercept = median.score[median.score$type == "non_eG4",2],linetype=2)
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"G4Structure.pdf"),
       device = "pdf",width = 6,height = 4)



#### GC含量 ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "001.2."
setwd("/home/yuss/flyG4/result/PQS")

path <- "/home/yuss/flyG4/result/PQS/GC-content"
files <- dir(path)
filespath <- lapply(files, function(x)paste(path,x,sep = '/'))
data <- list()
data <- lapply(filespath, function(x)fread(x))
a <- gsub(".GC.txt"," ",files) ##gsub("目标字符", "替换字符", 对象)
b <- gsub("001.2."," ",a)
data2 <- list()
for (i in 1:5) {
  data2[[i]] <- data[[i]][,c(1:6,15)]
  data2[[i]]$group <- b[i]
}
df <- do.call("rbind",data2)
#as.character(head(df$group)) 这列有空格
df$group <- gsub(" ","",df$group)
df$group %<>% factor(.,levels = c("non_eG4","kc_specific","s2_specific","overlap","merge")) ##设置顺序
colnames(df)[7] <- 'gc'
#*每类G4-------------------------------------------------------------------------------------
medians <- aggregate(gc ~  group, df, median)
my_comparisons = list(c("non_eG4","kc_specific"),c("kc_specific","s2_specific"),
                      c("s2_specific","overlap"),c("overlap","merge"))
b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c('#4d4d4d',b)
ggplot(data = df,aes(x = group,y = gc, color = group)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white",size = 0.8,fill = 'white') +
  stat_compare_means(comparisons = my_comparisons,label.y = c(rep(c(1,1.1),3)),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size = 4) +
  scale_color_manual(values = color ) +
  cowplot::theme_half_open() +
  theme(axis.text = element_text(size=14),
        axis.title.y = element_text(size = 16),
        axis.title.x =element_blank(),legend.position="none")+
  coord_cartesian(ylim = c(0.4,1.2)) +
  ylab("GC content") +
  geom_hline(yintercept = medians[medians$group =='non_eG4',2],linetype=2)
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"G4GCcontent.pdf"),
       device = "pdf",width = 6,height = 4)
#Ipaper::write_fig(p,file = "gc_content.group.pdf",width = 6,height = 4,devices = NULL,res = 300,show = F)

#*non_eG4和mergeG4的GC含量-----------------------------------------------------------------
gc.eG4.non_eG4 <- df[df$group=="non_eG4"|df$group=="merge"]
gc.eG4.non_eG4$group <- factor(gc.eG4.non_eG4$group,c("non_eG4","merge"))
my_comparisons = list(c("non_eG4","merge"))
ggplot(data = gc.eG4.non_eG4,aes(x = group,y = gc, fill = group)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,label.y = 0.95,
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size = 4) +
  scale_fill_manual(values = c("#A1A9D0","#F0988C") ) +
  cowplot::theme_half_open() +
  theme(axis.text = element_text(size=14),
        axis.title.y = element_text(size = 16),
        axis.title.x =element_blank(),legend.position="none")+
  coord_cartesian(ylim = c(0.4,1)) +
  ylab("GC content") +
  scale_x_discrete(labels = c("non-eG4","eG4"))
  # geom_hline(yintercept = medians[medians$group =='non_eG4',2],linetype=2)
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"non_eG4.eG4GCcontent.pdf"),
       device = "pdf",width = 3.2,height = 3.2)

#*overlap和细胞特异性G4的GC含量-----------------------------------------------------------------
gc.overlap.kc.s2 <- df[df$group=="overlap"|df$group=="kc_specific"|df$group=="s2_specific"]
gc.overlap.kc.s2$group <- factor(gc.overlap.kc.s2$group,levels = c("overlap","kc_specific","s2_specific"))
my_comparisons = list(c("overlap","kc_specific"),c("kc_specific","s2_specific"))
b = colorRampPalette(colors = c('#C18383', '#FDDBC7'))(3)
median.length <- aggregate(length ~  type, merge.all.3, median)
ggplot(data = gc.overlap.kc.s2,aes(x = group,y = gc, fill = group)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,label.y = c(0.92,0.95),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size = 4) +
  scale_fill_manual(values = b ) +
  cowplot::theme_half_open() +
  theme(axis.text = element_text(size=14),
        axis.title.y = element_text(size = 16),
        axis.title.x =element_blank(),legend.position="none")+
  coord_cartesian(ylim = c(0.4,1)) +
  ylab("GC content") +
  scale_x_discrete(labels = c("Overlap","Kc specific","S2 specific"))
# geom_hline(yintercept = medians[medians$group =='non_eG4',2],linetype=2)
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"OverlapAndSpecificG4.GCcontent.pdf"),
       device = "pdf",width = 4, height = 3.5)

rm(list = ls());gc();rm(list = ls())#清空
Num = "001.2."

#### 标准化后的density ####
#标准化后的density
merge.all <- fread("/home/yuss/flyG4/result/PQS/001.2.merge.all.bed") %>% as.data.frame()
chr.num <- table(merge.all$chr,merge.all$type)%>% data.frame()
colnames(chr.num) <- c("chr","type","Freq")
chr.num$chr <- factor(chr.num$chr,levels =c('2L','2R','3L','3R',4,'X','Y'))
chr.num$type <- factor(chr.num$type,levels = c("non_eG4","kc_specific","s2_specific","overlap","merge"))
chr.length <- fread("/home/yuss/flyG4/result/PQS/001.2.chr.length.txt") %>% data.frame()  
chr.gc <- fread("/home/yuss/flyG4/result/PQS/001.2.chr.GC.txt") %>% data.frame()
chr.num$gc <- chr.gc[match(chr.num$chr,chr.gc$X.1_usercol),5]
chr.num$length <- chr.length[match(chr.num$chr,chr.length$V1),2]
chr.num$density <- chr.num$Freq/(chr.num$gc*chr.num$length)
chr.num$density1000 <- chr.num$density*1000

b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c('#4d4d4d',b)
ggplot(chr.num, aes(x=chr,y=density1000,fill=type)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(),width = 0.9,color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值，#width设置矩形条的宽度
  coord_cartesian(ylim = c(0,0.89)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) + ##消除x轴与绘图区的间隙
  scale_fill_manual(values = color) +
  #geom_text(aes(label=Freq,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题
  ylab(bquote(Normalized~density~(x10^-3))) + 
  #ylab("Normalized density (No./kb/GC%)") +
  xlab("Chromosome arm") +
  #scale_fill_discrete(name="group") + ##修改图例标题名称
  theme(axis.title.y = element_text(size = 16), ##y坐标轴标题字体大小
        axis.title.x = element_text(size = 16), ##x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.position = c(0, 0.97), ##图例在绘图区域的位置 
        legend.direction = "horizontal", ##设置图例水平放置
        legend.title = element_blank()) ##删除图例标题名称
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"G4.ChrDensity.pdf"),
       device = "pdf",width = 6.5, height = 4.2)  

x.chr.num <- chr.num[chr.num$chr=="X",]
##分组求和,求每类G4在不同染色体的总数量
chr.num$type <- factor(chr.num$type,c("kc_specific","merge","non_eG4","overlap","s2_specific"))
x.chr.num$sum <- aggregate(chr.num$density1000, by=list(type=chr.num$type),sum)
##求每类G4在X染色体的比例：每类G4在X染色体的数量除以在不同染色体的总数量
x.chr.num$ratio <- x.chr.num$density1000/x.chr.num$sum$x
x.chr.num$percent <- round(x.chr.num$ratio*100, 2)
overlap.x.chr.num <- x.chr.num[c(1,4:5),]
overlap.x.chr.num$type <- factor(overlap.x.chr.num$type, levels = c("overlap","kc_specific","s2_specific"))

b = colorRampPalette(colors = c('#C18383', '#FDDBC7'))(3)
significance_labels <- data.frame(type = c("kc_specific", "s2_specific"),
  label = c("***", "*"))
ggplot(overlap.x.chr.num, aes(x=type,y=percent,fill=type)) +
  geom_bar(stat = 'identity',position = position_dodge(),width = 0.7,color = 'white')+
  coord_cartesian(ylim = c(0,50)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = b) +
  # geom_text(aes(label=percent,vjust=-0.5),color="black",size=4) +
  cowplot::theme_half_open() +
  ylab("(%)Normalized Proportion")+
  xlab("")+
  theme(axis.title.y = element_text(size = 16),
        axis.text = element_text(size=14),
        legend.position = "none")+
  scale_x_discrete(labels = c("Overlap","Kc specific","S2 specific")) +
  # geom_signif(y_position=30, xmin=2, xmax=4,
  #             annotation=c("***"), tip_length=0.06)#p-value < 2.2e-16
  geom_text(
    data = significance_labels,
    aes(x = type, y = c(40,11), label = label),
    vjust = -0.5,
    size = 5
  )
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"Normalized.OverlapAndSpecificG4.XchrProportion.pdf"),
       device = "pdf",width = 4.4, height = 4.4)

##显著性
x.chr.num$density1000.chA <- x.chr.num$sum$x - x.chr.num$density1000
x.chr.num$type
x.chr.num$density1000
x.chr.num$density1000.chA

data1 <- matrix(c(1582,2452,975,7595), nrow = 2)
colnames(data1) <- c("kc","overlap")
rownames(data1) <- c("chr X","chr A")
fisher.test(data1) ##p-value < 2.2e-16

data2 <- matrix(c(231,2115,975,7595), nrow = 2)
colnames(data2) <- c("s2","overlap")
rownames(data2) <- c("chr X","chr A")
fisher.test(data2) ##p-value = 0.03737

data1 <- matrix(c(1582,2452,1206,9710), nrow = 2)
colnames(data1) <- c("kc","other")
rownames(data1) <- c("chr X","chr A")
fisher.test(data1) ##p-value < 2.2e-16

data2 <- matrix(c(231,2115,2557,10047), nrow = 2)
colnames(data2) <- c("s2","other")
rownames(data2) <- c("chr X","chr A")
fisher.test(data2) ##p-value < 2.2e-16

#### 保守性 ####
#*每类G4的保守性--------------------------------------------------------------------
merge.all <- fread("/home/yuss/flyG4/result/PQS/001.2.merge.all.bed") %>% data.frame()
phastCons <- fread("/home/yuss/flyG4/result/PQS/001.2.phastCons.bed") %>% data.frame()
phyloP <- fread("/home/yuss/flyG4/result/PQS/001.2.phyloP.bed") %>% data.frame()  
merge.all$pha <- phastCons[match(merge.all$id,phastCons$V1),6]
merge.all$phy <- phyloP[match(merge.all$id,phyloP$V1),6]
merge.all$type <- factor(merge.all$type,levels =c("non_eG4","kc_specific","s2_specific","overlap","merge"))
my_comparisons = list(c("non_eG4","kc_specific"),c("kc_specific","s2_specific"),
                      c("s2_specific","overlap"),c("overlap","merge"))
b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c('#4d4d4d',b)
median.pha <- aggregate(pha ~  type, merge.all, median)
library(ggpubr)
ggplot(data = merge.all,aes(x=type,y=pha,fill=type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,label.y = c(rep(c(1,1.1),3)),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = color) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("PhastCons score") +
  #coord_cartesian(ylim = c(0,80)) +
  geom_hline(yintercept = median.pha[median.pha$type == "non_eG4",2],linetype=2)
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"G4PhastCons.pdf"),
       device = "pdf",width = 6,height = 4)
median.phy <- aggregate(phy ~  type, merge.all, median)
ggplot(data = merge.all,aes(x=type,y=phy,fill=type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,label.y = c(rep(c(1.9,2.1),3)),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = color) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("PhyloP score") +
  coord_cartesian(ylim = c(0,2.5)) +
  geom_hline(yintercept = median.phy[median.phy$type == "non_eG4",2],linetype=2)
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"G4PhyloP.pdf"),
       device = "pdf",width = 6, height = 4)  

#*non_eG4和eG4--------------------------------------------------------------------------------------
non_eG4.eG4 <- merge.all[merge.all$type=="non_eG4"|merge.all$type=="merge",]
my_comparisons = list(c("non_eG4","merge"))
ggplot(data = non_eG4.eG4,aes(x=type,y=pha,fill=type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,label.y = 1,
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#A1A9D0","#F0988C")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("PhastCons score") +
  scale_x_discrete(labels = c("non-eG4", "eG4"))
#geom_hline(yintercept = median.length[median.length$type == "non_eG4",2],linetype=2)
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"non_eG4.eG4PhastCons.pdf"),
       device = "pdf",width = 3.5, height = 3.5) 

ggplot(data = non_eG4.eG4,aes(x=type,y=phy,fill=type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,label.y = 1.6,
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#A1A9D0","#F0988C")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("PhyloP score") +
  coord_cartesian(ylim = c(0,2.0)) +
  scale_x_discrete(labels = c("non-eG4", "eG4"))
#geom_hline(yintercept = median.length[median.length$type == "non_eG4",2],linetype=2)
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"non_eG4.eG4PhyloP.pdf"),
       device = "pdf",width = 3.5, height = 3.5) 

#*overlap和细胞特异性G4的保守性--------------------------------------------------------------------------------------
overlap <- merge.all[merge.all$type=="overlap"|merge.all$type=="kc_specific"|merge.all$type=="s2_specific",]
my_comparisons = list(c("overlap","kc_specific"),c("kc_specific","s2_specific"))
overlap$type <- factor(overlap$type,c("overlap","kc_specific","s2_specific"))

b = colorRampPalette(colors = c('#C18383', '#FDDBC7'))(3)
ggplot(data = overlap,aes(x=type,y=pha,fill=type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,label.y = c(1.02,1.1),
                     tip.length = 0,size=4)+
  scale_fill_manual(values = b) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("PhastCons score") +
  coord_cartesian(ylim = c(0,1.2)) +
  scale_x_discrete(labels = c("Overlap","Kc specific","S2 specific"))
#geom_hline(yintercept = median.length[median.length$type == "overlap",2],linetype=2)
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"OverlapAndSpecificG4.PhastCons.pdf"),
       device = "pdf",width = 4, height = 3.5)

ggplot(data = overlap,aes(x=type,y=phy,fill=type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,label.y = c(1.82,1.90),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = b) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("PhyloP score") +
  coord_cartesian(ylim = c(0,2.2)) +
  scale_x_discrete(labels = c("Overlap","Kc specific","S2 specific"))
#geom_hline(yintercept = median.length[median.length$type == "overlap",2],linetype=2)
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"OverlapAndSpecificG4.PhyloP.pdf"),
       device = "pdf",width = 4, height = 3.5)

#### 同源 ####
#人与果蝇同源的G4占人中G4的比例
BiocManager::install("liftOver")

#### ####
eG4 <- fread("/home/yuss/flyG4/result/PQS/001.2.merge.bed") %>% as.data.frame()
eG4$length <- eG4$end-eG4$start
sum(eG4$length)
#502023
#基因组总长度143,726,002