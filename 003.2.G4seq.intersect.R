rm(list = ls());gc();rm(list = ls())#清空
Num = "003.2."
##读取数据
pqs.k <- fread("/home/yuss/flyG4/result/data.reliability/003.2.pqs.G4-seqK.bed") %>% as.data.frame()
pqs.pds <- fread("/home/yuss/flyG4/result/data.reliability/003.2.pqs.G4-seqPDS.bed") %>% as.data.frame()
pqs.G4seq <- bind_cols(pqs.k,pqs.pds$V7)
colnames(pqs.G4seq) <- c("chr","start","end","pqs.id","score","strand","k","pds")
pqs.G4seq$k <- ifelse(pqs.G4seq$k==0,0,1)
pqs.G4seq$pds <- ifelse(pqs.G4seq$pds==0,0,2)
pqs.G4seq$sum <- rowSums(pqs.G4seq[,7:8])
table(pqs.G4seq$sum)
no_G4_formed <- pqs.G4seq[pqs.G4seq$sum==0,]
k <- pqs.G4seq[pqs.G4seq$sum==1,]
pds <- pqs.G4seq[pqs.G4seq$sum==2,]
overlap <- pqs.G4seq[pqs.G4seq$sum==3,]
merge <- pqs.G4seq[pqs.G4seq$sum!=0,]
no_G4_formed$type <- "No_G4_formed"
k$type <- "K+"
pds$type <- "PDS"
overlap$type <- "overlap"
merge$type <- "merge"
merge.all <- bind_rows(no_G4_formed,k,pds,overlap,merge) 
merge.all$type <- factor(merge.all$type,levels = c("No_G4_formed","K+","PDS","overlap","merge"))

##保存
write.table(k,file = paste0("/home/yuss/flyG4/result/data.reliability/",Num,"k+.bed"),sep = '\t',
            col.names = F,row.names = F,quote = F)
write.table(pds,file = paste0("/home/yuss/flyG4/result/data.reliability/",Num,"PDS.bed"),sep = '\t',
            col.names = F,row.names = F,quote = F)
write.table(overlap, file = paste0("/home/yuss/flyG4/result/data.reliability/",Num,"overlap.bed"),sep = '\t',
            col.names = F,row.names = F)
write.table(merge, file = paste0("/home/yuss/flyG4/result/data.reliability/",Num,"merge.bed"),sep = '\t',
            col.names = F,row.names = F)
write.table(no_G4_formed, file = paste0("/home/yuss/flyG4/result/data.reliability/",Num,"no_G4_formed.bed"),sep = '\t',
            col.names = F,row.names = F,quote = F)
#### 结构稳定性 ####
median.score <- aggregate(score ~ type, merge.all, median)
my_comparisons = list(c("No_G4_formed","K+"),c("K+","PDS"),
                      c("PDS","overlap"),c("overlap","merge"))
b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c('#4d4d4d',b)
library(ggpubr)
ggplot(data = merge.all,aes(x=type,y=score,fill=type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,label.y = c(rep(c(172,185),3)),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4) +
  scale_fill_manual(values = color) +
  cowplot::theme_half_open() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Structural stability") +
  coord_cartesian(ylim = c(50,200)) +
  geom_hline(yintercept = median.score[median.score$type == "No_G4_formed",2],linetype=2)

tapply(merge.all$score, merge.all$type, median)
table(merge.all$type)
median(merge.all$score)

#### G4基因组注释 ####
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
Dmel.txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
#peaks注释,tss(转录起始位点)选取上下游500bp，可自行调整(参考文章中的取值)
anno <- annotatePeak("/home/yuss/flyG4/result/data.reliability/003.2.no_G4_formed.bed",
                     tssRegion = c(-500,500),TxDb = Dmel.txdb)
anno2 <- annotatePeak("/home/yuss/flyG4/result/data.reliability/003.2.k+.bed",
                      tssRegion = c(-500,500),TxDb = Dmel.txdb)
anno3 <- annotatePeak("/home/yuss/flyG4/result/data.reliability/003.2.PDS.bed",
                      tssRegion = c(-500,500),TxDb = Dmel.txdb)
anno4 <- annotatePeak("/home/yuss/flyG4/result/data.reliability/003.2.overlap.bed",
                      tssRegion = c(-500,500),TxDb = Dmel.txdb)
anno5 <- annotatePeak("/home/yuss/flyG4/result/data.reliability/003.2.merge.bed",
                      tssRegion = c(-500,500),TxDb = Dmel.txdb)
df1 <- anno@annoStat %>% as.data.frame()
df2 <- anno2@annoStat %>% as.data.frame()
df3 <- anno3@annoStat %>% as.data.frame()
df4 <- anno4@annoStat %>% as.data.frame()
df5 <- anno4@annoStat %>% as.data.frame()
df1$Type <- 'no_G4_formed'
df2$Type <- 'k+'
df3$Type <- 'pds'
df4$Type <- 'overlap'
df5$Type <- 'merge'
df <- rbind(df1,df2,df3,df4,df5)
df$Type %<>% factor(.,levels = c("no_G4_formed","k+","pds","overlap","merge"))

b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(5)
color <- c('#4d4d4d',b)
#画图
library(ggplot2)
ggplot(df,aes(x = Feature,y= Frequency,fill = Type)) +
  geom_bar(position = 'dodge',stat = 'identity',color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值
  scale_fill_manual(values = color) +
  scale_y_continuous(expand = c(0,0)) +
  cowplot::theme_half_open() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_blank(),
        axis.text.x = element_text(vjust = 0.5,hjust = 1,size = 14), #hjust = 1右对齐
        axis.text.y = element_text(size = 14),
        axis.ticks.x = element_blank(),
        legend.position = 'top') +
  ylab("Frequency (%)") +
  geom_text(aes(label=str_c(round(Frequency,2),"%"),y=Frequency + 3), color="black", size=3.5,
                position = position_dodge(.9))+ #y=Frequency + 3是高度，str_c函数是将多个字符串合并成一个字符串，y=Frequency + 3是标签高低位置设置，position = position_dodge(.9)是标签左右位置设置
  coord_flip(ylim = c(0,40)) #x 和 y 翻转的笛卡尔坐标

rm(list = ls());gc();rm(list = ls())#清空
Num = "003.2."
##读取数据
pqs.G4seqmerge <- fread("/home/yuss/flyG4/result/data.reliability/003.2.pqs.G4seqmerge.bed") %>% as.data.frame()
colnames(pqs.G4seqmerge) <- c("chr","start","end","pqs_id","score","strand","k+PDS")
pqs.G4seqmerge$type <- ifelse(pqs.G4seqmerge$`k+PDS`==0,"No_G4_formed","G4")
pqs.G4seqmerge$type <- factor(pqs.G4seqmerge$type,levels=c("No_G4_formed","G4"))
median.score <- aggregate(score ~ type, pqs.G4seqmerge, median)
pqs.G4seqmerge.No_G4_formed <- pqs.G4seqmerge[pqs.G4seqmerge$`k+PDS`==0,]
pqs.G4seqmerge.G4 <- pqs.G4seqmerge[pqs.G4seqmerge$`k+PDS`!=0,]
write.table(pqs.G4seqmerge.No_G4_formed,file = paste0("/home/yuss/flyG4/result/data.reliability/",Num,"pqs.G4seqmerge.No_G4_formed"),sep = '\t',
            col.names = F,row.names = F,quote = F)
write.table(pqs.G4seqmerge.G4,file = paste0("/home/yuss/flyG4/result/data.reliability/",Num,"pqs.G4seqmerge.G4.bed"),sep = '\t',
            col.names = F,row.names = F,quote = F)

#### 结构稳定性 ####
my_comparisons = list(c("No_G4_formed","G4"))
library(ggpubr)
ggplot(data = pqs.G4seqmerge,aes(x=type,y=score,fill=type))+
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,label.y = 130,
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4) +
  scale_fill_manual(values = c("#4d4d4d",'#C18383')) +
  cowplot::theme_half_open() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Structural stability") +
  coord_cartesian(ylim = c(50,150)) +
  scale_x_discrete(labels = c("No_G4_formed" = "non_G4", "G4" = "G4")) +
  geom_hline(yintercept = median.score[median.score$type == "No_G4_formed",2],linetype=2)  
ggsave(filename = paste0("/home/yuss/flyG4/result/data.reliability/Picture/",Num,"G4seq.StructuralStability.pdf"),
       device = "pdf",width = 6,height = 4)

#### G4基因组注释 ####
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
Dmel.txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
annoG4 <- annotatePeak("/home/yuss/flyG4/result/data.reliability/003.2.pqs.G4seqmerge.G4.bed",
                       tssRegion = c(-500,500),TxDb = Dmel.txdb)
annoNO_G4_formed <- annotatePeak("/home/yuss/flyG4/result/data.reliability/003.2.pqs.G4seqmerge.No_G4_formed",
                                 tssRegion = c(-500,500),TxDb = Dmel.txdb)
df1 <- annoG4@annoStat %>% as.data.frame()
df2 <- annoNO_G4_formed@annoStat %>% as.data.frame()
df1$Type <- 'G4'
df2$Type <- 'NO_G4_formed'
df <- rbind(df1,df2)
df$Type %<>% factor(.,levels = c("NO_G4_formed","G4"))
library(ggplot2)
ggplot(df,aes(x = Feature,y= Frequency,fill = Type)) +
  geom_bar(position = 'dodge',stat = 'identity',color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值
  scale_fill_manual(values = c("#4d4d4d",'#C18383')) +
  scale_y_continuous(expand = c(0,0)) +
  cowplot::theme_half_open() +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_blank(),
        axis.text.x = element_text(vjust = 0.5,hjust = 1,size = 14), #hjust = 1右对齐
        axis.text.y = element_text(size = 14),
        axis.ticks.x = element_blank(),
        legend.position = 'top') +
  ylab("Frequency (%)") +
  geom_text(aes(label=str_c(round(Frequency,2),"%"),y=Frequency + 3), color="black", size=3.5,
            position = position_dodge(.9))+ #y=Frequency + 3是高度，str_c函数是将多个字符串合并成一个字符串，y=Frequency + 3是标签高低位置设置，position = position_dodge(.9)是标签左右位置设置
  coord_flip(ylim = c(0,40)) #x 和 y 翻转的笛卡尔坐标
