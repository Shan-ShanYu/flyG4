#### 图3-13 eG4 结构按组织特异性分组####
##合并数据，整理表格，然后对sum进行table,然后根据G4在几种细胞系中存在进行分类，分成6个group,并画图，设置x = Var1,y= num，把颜色映射给分组，再利用coord_flip翻转xy轴
# 人两次 G4_pqs数据合并-------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())
library(magrittr)
library(data.table)
library(tidyverse)

df1 <- fread('/NAS/yulix/all_G4_data/01.G4_pqs/hg_G4_cell.bed') %>% data.frame()
df2 <- fread('/NAS/yulix/all_G4_data/01.G4_pqs/hg.pqs.result') %>% data.frame()

df1 <- df1[,c(4,7:13)]
df2 <- df2[,-65]
df3 <- left_join(df2,df1,by = c('id' = 'id'))

colnames(df3)[c(8,65)]
colnames(df3)[c(31,67)]
colnames(df3)[c(23,64,71)]

df4 <- df3[,-c(8,65,31,67,23,64,71)]

a549 <- df3[,c(4,8,65)]
HepG2 <- df3[c(4,31,67)]
HEK293T <- df3[c(4,23,64,71)]

df4$A549 <- rowSums(a549[,2:3])
df4$A549 <- if_else(df4$A549 >= 1,1,0)

df4$HepG2 <- rowSums(HepG2[,2:3])
df4$HepG2 <- if_else(df4$HepG2 >= 1,1,0)

df4$HEK293T <- rowSums(HEK293T[,2:4])
df4$HEK293T <- if_else(df4$HEK293T >= 1,1,0)

df4$sum <- rowSums(df4[,7:67])

write.table(df4,file = '/NAS/yulix/all_G4_data/01.G4_pqs/hg_G4_pqs.txt',
            sep = '\t',row.names = F,quote = F)

df4 <- fread('/NAS/yulix/all_G4_data/01.G4_pqs/hg_G4_pqs.txt') %>% data.frame()
pqs <- df4[df4$sum == 0,1:6]
G4 <- df4[df4$sum != 0,1:6]
write.table(G4,file = '/NAS/yulix/all_G4_data/01.G4_pqs/hg_G4.bed',
            sep = '\t',row.names = F,col.names = F,quote = F)
write.table(pqs,file = '/NAS/yulix/all_G4_data/01.G4_pqs/hg_pqs.bed',
            sep = '\t',row.names = F,col.names = F,quote = F)
# 人 G4 数量
sum(df4$sum != 0)

# number_of_G4-------------------------------------------------------------------------
a <- table(df4$sum) %>% data.frame()
a <- a[-1,]
colnames(a) <- c('Var1','num')

a$Var1= factor(a$Var1,levels = factor(order(-(as.numeric(a$Var1))))) ##order()函数对数据框进行排序,返回的值表示位置，默认是升序，依次对应的是向量的最小值、次小值、第三小值…最大值
a$group <-if_else(a$Var1 ==1,'Group1','NA')
a$group <-if_else(a$Var1 ==2,'Group2',a$group)
a$group <-if_else(a$Var1 ==3,'Group3',a$group)
a$group <-if_else(a$Var1 == 4 | a$Var1 == 5,'Group4',a$group)
a$group <-if_else(a$Var1 == 6 | a$Var1 == 7 | a$Var1 == 8
                  | a$Var1 == 9 | a$Var1 == 10,'Group5',a$group)
a$group <-if_else(a$group =='NA','Group6',a$group)

b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)

library(ggplot2)
p <- ggplot(data = a,aes(x = Var1,y= num)) +
  geom_col(aes(fill = group),position=position_dodge(1),width = 1,color = 'white') +
  coord_flip(ylim = c(0,140000)) + #
  scale_y_continuous(expand = c(0,0)) +
  geom_text(aes(label=num,y=num + 3500), color="black", size=5) +
  cowplot::theme_half_open()+
  scale_fill_manual(values = b) +
  ylab("Number of eG4s ") + 
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=16),
        axis.ticks.x = element_blank(),
        legend.position = "none")

library(Ipaper,lib.loc = "/home/yulix/R/x86_64-pc-linux-gnu-library/4.0/")
write_fig(p,file = "/NAS/yulix/all_G4_data/figure/number_of_G4.pdf",
          width = 8,height = 12,devices = NULL,res = 300,show = F) 

#### 图3-14 G4 CUT&Tag样本的eG4 结构的组织特异性形成 ####
##人G4 CUT&Tag数据进行table，再x = Var1,y= num,fill = Var1画图
df <- fread("/NAS/yulix/20221102cut-tag/result/peak_deal/human/hg_G4_cell.bed") %>% data.frame()
b <- table(df$sum) %>% data.frame()
colnames(b) <- c('Var1','num')
b <- b[-1,]

b$Var1= factor(b$Var1,levels = factor(order(-(as.numeric(b$Var1)))))

p <- ggplot(data = b,aes(x = Var1,y= num,fill = Var1)) +
  geom_bar(stat = 'identity',position=position_dodge(0.7),width = 0.7,
           color = 'white',fill = "#e1b1bd") +
  scale_y_continuous(expand = c(0,0)) +
  geom_text(aes(label=num),hjust = 1, color="black", size=4) +
  cowplot::theme_half_open()+
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_blank(),
        axis.text = element_text(size=14),
        legend.position = "none")+
  ylab("Number of eG4s in Human")+
  coord_flip(ylim = c(0,42000))

library(Ipaper,lib.loc = "/home/yulix/R/x86_64-pc-linux-gnu-library/4.0/")
write_fig(p,file = "/NAS/yulix/all_G4_data/figure/Number_of_eG4s_in_Human.pdf",
          width = 5,height = 4,devices = NULL,res = 300,show = F) 

####图3-15 G4 CUT&Tag原始数据和G4 ChIP-seq公共数据之间eG4 组织特异性的相关性####
##先计算
df1 <- fread('/NAS/yulix/all_G4_data/01.G4_pqs/hg_G4_cell.bed') %>% data.frame()
df2 <- fread('/NAS/yulix/all_G4_data/01.G4_pqs/hg_G4_pqs.txt') %>% data.frame()

df2 <- df2[df2$sum != 0,]
a <- df2[,1:4]
a$x <- df1[match(a$id,df1$id),14]
a$y <- df2[match(a$id,df2$id),68]
a$x2 <- (7-a$x)/7  ##推测计算的是没有交集的比例
a$y2 <- (53-a$y)/53

cor.test(a$x2,a$y2)
library(ggpointdensity,lib.loc = "/home/yulix/R/x86_64-pc-linux-gnu-library/4.0/")
library(Ipaper,lib.loc = "/home/yulix/R/x86_64-pc-linux-gnu-library/4.0/")

P <- ggplot(data = a,aes(x = x2,y = y2)) +
  #geom_pointdensity(shape = 21,size = 1,color = '#486c8d',fill = '#486c8d') +
  #geom_point(shape = 21,size = 1,color = '#486c8d',fill = '#486c8d') +
  theme_bw()+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 14))+
  xlab(paste0("Cell type specificity of eG4s","\nin our G4 CUT&Tag data"))+ #\n换行
  ylab(paste0("Cell type specificity of eG4s","\nin public G4 ChIP-seq data"))+
  annotate("text",x = 0.2,y = 0.95,label = "R = 0.76; P < 2e-16",size = 5) + #在(0.2,0.95)处添加标签
  stat_smooth(method = 'gam',formula = y~x, level = 0.95,color = "#b2182b") #stat_smooth绘制线性拟合线,level = 0.95表示95%的置信区间,不添加置信区间se = FALSE

Ipaper::write_fig(p,file = "/NAS/yulix/all_G4_data/figure/TissueSpecificity.pdf",
                  width = 6,height = 4.5,devices = NULL,res = 300,show = F) 

a2 <- a[,7:8]
a2 <- unique(a2)
p <- ggplot(data = a2,aes(x = x2,y = y2)) +
  geom_point(shape = 21,size = 1,color = '#486c8d',fill = '#486c8d') +
  theme_bw()+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 14))+
  xlab("Tissue specificity of eG4s in our G4 CUT&Tag data")+
  ylab("Tissue specificity of eG4s in public G4 ChIP-seq data")+
  annotate("text",x = 0.2,y = 0.95,label = "R = 0.76; P < 2e-16",size = 5)

Ipaper::write_fig(p,file = "/NAS/yulix/all_G4_data/figure/TissueSpecificity2.pdf",
                  width = 6,height = 4.5,devices = NULL,res = 300,show = F) 

##绘制散点图时，若散点数目很多，散点之间相互重叠，则不易观察散点趋势，此时可绘制密度散点图解决
library(ggpointdensity,lib.loc = "/home/yulix/R/x86_64-pc-linux-gnu-library/4.0/")

b <- head(a,1000)
p <- ggplot(data = a,aes(x = x2,y = y2)) +
  geom_pointdensity() +
  #geom_point(shape = 21,size = 1,color = '#486c8d',fill = '#486c8d') +
  theme_bw()+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),legend.position = "none")+
  xlab(paste0("Cell type specificity of eG4s","\nin our G4 CUT&Tag data"))+
  ylab(paste0("Cell type specificity of eG4s","\nin public G4 ChIP-seq data"))+
  annotate("text",x = 0.9,y = 0.4,label = paste0("R = 0.76","\nP < 2e-16"),size = 5) +
  stat_smooth(method = 'gam',formula = y~x, se = TRUE, level = 0.95,color = "#b2182b")+
  scale_color_gradient(low = "grey", high = "brown") ##设置双色梯度

Ipaper::write_fig(p,file = "/NAS/yulix/all_G4_data/figure/TissueSpecificity.pdf",
                  width = 5,height = 4,devices = NULL,res = 300,show = F) 

p <- ggplot(data = a,aes(x = x2,y = y2)) +
  geom_pointdensity()+
  theme_bw()+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),legend.position = "none")+
  xlab(paste0("Cell type specificity of eG4s","\nin our G4 CUT&Tag data"))+
  ylab(paste0("Cell type specificity of eG4s","\nin public G4 ChIP-seq data"))+
  annotate("text",x = 0.2,y = 0.95,label = "R = 0.76; P < 2e-16",size = 5)+
  stat_smooth(method = 'gam',formula = y~x, se = TRUE, level = 0.95,color = "black")+
  scale_color_gradient(low = "grey", high = "brown")

ggpointdensity::stat_density2d(geom = "polygon")
geom_bin2d()

####图3-17 eG4 结构的组织特异性形成####
# 小鼠的peak 与pqs交集-------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())
library(magrittr)
library(data.table)
library(tidyverse)

path <- '/NAS/yulix/all_G4_data/01.G4_pqs/mm.intersect'
files <- dir(path)
filespath <- lapply(files, function(x)paste(path,x,sep = '/'))
data <- list()
data <- lapply(filespath,function(x)fread(x))

a = files
a %<>% gsub('1001-','',.) %>% gsub('_R1_val_1','',.) %>% gsub('_peaks.txt','',.)%>% gsub('0918-','',.) ##把多余的名字符号删除

data2 <- data.frame(matrix(NA,nrow(data[[1]]),(length(data)+6))) #设置一个行为1516050，列为33的空表格
data2[,1:6] = data[[1]][,1:6] #前六列的信息
for (i in 1:length(data)){
  data2[,i+6] = data[[i]]$V7
} #把所有细胞的交集信息以循环的形式附在表格后边
colnames(data2) <- c('chr','start','end','id','score','strand',a)

df = data2
b1 <- lapply(colnames(df)[c(7:12)], function(x)strsplit(x,'-',fixed = T)[[1]][1]) %>% unlist
b2 <- lapply(colnames(df)[c(13:17)], function(x)strsplit(x,'_',fixed = T)[[1]][1]) %>% unlist
b3 <- colnames(df)[c(18:33)] %>% gsub('_CUT-Tag','',.)
b3 <- c("mESC_DMSO",          "mESC_DMSO" ,        
        "mESC_DRB",           "mESC_DRB" ,         
        "mESC_mock",          "mESC_Mungbean-100U",
        "mESC_Mungbean-150U", "mESC",              
        "mESC",               "mESC",              
        "mESC",               "mESC",              
        "mESC_triptolide",    "mESC_triptolide",   
        "mNPC",               "mNPC")
b <- c(b1,b2,b3)

colnames(df)[7:ncol(df)] <- b
rownames(df) <- df$id
c = apply(df[,7:ncol(df)], 1,function(x)tapply(x,b,sum))
d = data.frame(t(as.matrix(c)))

d$mESC_mock <- d$mESC_mock +1
d$`mESC_Mungbean.100U` <- d$`mESC_Mungbean.100U` + 1 
d$`mESC_Mungbean.150U`  <- d$`mESC_Mungbean.150U` +1

for (i in 1:ncol(d)){
  d[,i] <- if_else(d[,i] >= 2,1,0)
}

d$sum <- rowSums(d)
df2 <- data.frame(df[,1:6],d)
write.table(df2,
            file="/NAS/yulix/all_G4_data/01.G4_pqs/mm_G4_pqs.txt",
            sep = '\t',row.names = FALSE,col.names = T,quote=FALSE)

g4 <- df2[df2$sum !=0,1:6]
write.table(g4,
            file="/NAS/yulix/all_G4_data/01.G4_pqs/mm_G4.bed",
            sep = '\t',row.names = FALSE,col.names = F,quote=FALSE)
pqs <- df2[df2$sum ==0,1:6]
write.table(pqs,
            file="/NAS/yulix/all_G4_data/01.G4_pqs/mm_pqs.bed",
            sep = '\t',row.names = FALSE,col.names = F,quote=FALSE)

# Number of eG4s in Mouse-------------------------------------------------------------------------
a <- table(df2$sum) %>% data.frame()
a <- a[-1,]
colnames(a) <- c('Var1','num')

a$Var1= factor(a$Var1,levels = factor(order(-(as.numeric(a$Var1)))))

library(ggplot2)
p <- ggplot(data = a,aes(x = Var1,y= num)) +
  geom_col(fill = "#a5b586",position=position_dodge(1),width = 1,color = 'white') +
  coord_flip(ylim = c(0,22000)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_text(aes(label=num,y=num + 1000), color="black", size=4) +
  cowplot::theme_half_open()+
  ylab("Number of eG4s in Mouse") + 
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=16),
        axis.ticks.x = element_blank(),
        legend.position = "none")

library(Ipaper,lib.loc = "/home/yulix/R/x86_64-pc-linux-gnu-library/4.0/")
write_fig(p,file = "/NAS/yulix/all_G4_data/figure/Number_of_eG4s_in_Mouse.pdf",
          width = 6,height = 4,devices = NULL,res = 300,show = F)    

# Number of eG4s in Pig-------------------------------------------------------------------------
df <- fread('/NAS/yulix/all_G4_data/01.G4_pqs/pig_G4_pqs.txt') %>% data.frame()

b <- table(df$sum) %>% data.frame()
colnames(b) <- c('Var1','num')
b <- b[-1,]

b$Var1= factor(b$Var1,levels = factor(order(-(as.numeric(b$Var1))))) 

p <- ggplot(data = b,aes(x = Var1,y= num,fill = Var1)) +
  geom_bar(stat = 'identity',position=position_dodge(0.7),width = 0.7,
           color = 'white',fill = "#e1b1bd") +
  scale_y_continuous(expand = c(0,0)) +
  geom_text(aes(label=num),hjust = 1, color="black", size=5) +
  cowplot::theme_half_open()+
  theme(axis.title.x = element_text(size=14),
        axis.title.y = element_blank(),
        axis.text = element_text(size=12),
        legend.position = "none")+
  ylab("Number of eG4s in Pig")+
  coord_flip(ylim = c(0,30000))

library(Ipaper,lib.loc = "/home/yulix/R/x86_64-pc-linux-gnu-library/4.0/")
write_fig(p,file = "/NAS/yulix/all_G4_data/figure/Number_of_eG4s_in_pig.pdf",
          width = 4,height = 3,devices = NULL,res = 300,show = F)

#### 图3-23 每组中结构保守的eG4s的比例 ####
##计算hm（人、小鼠同源的eG4s），hp（人、猪同源的eG4s），hpm（人、小鼠、猪同源的eG4s）的数量，再除以人源eG4s的数量
rm(list = ls());gc();rm(list = ls())
hp <- fread("/NAS/yulix/all_G4_data/01.G4_pqs/human-pig") %>% data.frame()
hm <- fread("/NAS/yulix/all_G4_data/01.G4_pqs/human-mouse") %>% data.frame()
hpm <- fread("/NAS/yulix/all_G4_data/01.G4_pqs/human-pig-mouse") %>% data.frame()

list <- list(hp = hp,hm = hm,hpm = hpm)
for (i in 1:length(list)){
  list[[i]] <- list[[i]][list[[i]]$V8 != 0,]
  list[[i]]$class <- names(list)[i]
} 

df <- do.call('rbind',list) ##合并list,把表竖起来
a <- aggregate(df$V7,by = list(df$class),table)  %>% data.frame() ##按照df$class分组,分成了三组，对每组的group进行计总数

a <- tapply(df$V7,df$class,table) 
for (i in 1:length(a)){
  a[[i]] <- t(a[[i]]) %>% data.frame()
  a[[i]]$class <- names(a)[i]
}

b <- do.call('rbind',a)
group = c(paste0('Group',1:6))
num = c(128233,
        75258,
        51615,
        59420,
        52888,
        37739)# 6组人源eG4s的数量
df2 <- data.frame(group = group,num = num)
b$num <- df2[match(b$Var2,df2$group),2]
b$ratio <- b$Freq/b$num

p <- ggplot(b ,aes(x = Var2,y = ratio,fill = class)) +
  geom_bar(stat="identity",width= 0.75,position= position_dodge(0.75),color = 'white') +
  scale_fill_manual(values = c('#a5b586','#86a5b5','#e1b1bd')) +
  theme_bw()+
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y =element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 16),
        legend.position = c(0.1,0.85),legend.title = element_blank()) +
  ylab(paste0('(%) Conserved at structure level'))+
  coord_cartesian(ylim = c(0,0.35))

Ipaper::write_fig(p,file = "/NAS/yulix/all_G4_data/figure/Conserved.pdf",
                  width = 6,height = 5,devices = NULL,res = 300,show = F) 

#### 图3-24 人和小鼠保守的eG4s组织特异性的相关性 ####
##首先找到人和小鼠保守的eG4s是哪些，再找到交集信息，最后计算相关性
rm(list = ls());gc();rm(list = ls())
df1 <- fread('/NAS/yulix/all_G4_data/01.G4_pqs/mengwei.mm10-hg.G4.txt') %>% data.frame()
df2 <- fread('/NAS/yulix/all_G4_data/01.G4_pqs/mm_G4_pqs.txt') %>% data.frame()
df3 <- fread('/NAS/yulix/all_G4_data/01.G4_pqs/hg_G4_pqs.txt',nThread = 10) %>% data.frame()

a <- df1[,c(4,10)]
a$x <- df2[match(a$mm_id,df2$id),19] #分别把小鼠和人的有无交集的信息匹配过来
a$y <- df3[match(a$hg19_id,df3$id),68]

a$x2 <- (12-a$x)/12
a$y2 <- (52-a$y)/52

cor.test(a$x2,a$y2)

cor = cor.test(a$x2,a$y2)$estimate #相关性分数

p <- ggplot(data = a,aes(x = x2,y = y2)) +
  geom_point(shape = 21,size = 1,color = '#486c8d',fill = '#486c8d') +
  theme_bw()+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14))+
  xlab(paste0("Tissue specificity of orthologous eG4s","\nin mouse"))+
  ylab(paste0("Tissue specificity of orthologous eG4s","\nin human"))+
  annotate("text",x = 0.2,y = 0.9,label = paste0("R = ",round(cor,2),"; P < 2e-16"),size = 5) +
  stat_smooth(method = 'gam',formula = y~x, se = TRUE, level = 0.95,color = "#b2182b")

Ipaper::write_fig(p,file = "/NAS/yulix/all_G4_data/figure/tissueSpecificitymouse.pdf",
                  width = 6,height = 4.5,devices = NULL,res = 300,show = F) 

#### 图3-25 non-eG4s和eG4s的基因组分布B ####
##在启动子区(TSS)附近的 ChIP peaks 的可视化
rm(list = ls());gc();rm(list = ls())
path <- "/NAS/yulix/all_G4_data/01.G4_pqs/group-bed"
files <- dir(path)
filepath <- sapply(files, function(x){
  paste(path,x,sep = '/')
})

#加载包
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)

peak=list(`Non-eG4` = filepath[[7]],Group1= filepath[[1]],Group2 = filepath[[2]],
          Group3 = filepath[[3]],Group4 = filepath[[4]],Group5 = filepath[[5]],
          Group6 = filepath[[6]])

promoter <- getPromoters(TxDb=txdb, upstream=1500, downstream=1500) ##getPromoters函数定义启动子区域的范围
tagMatrixList <- lapply(peak, getTagMatrix, windows=promoter) ##getTagMatrix函数将Peaks匹配到启动子区域，然后计算出tagMatrix用于画图
tagHeatmap(tagMatrixList, xlim = c(-1500,1500),color = "red") ##热图
b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(6)
color <- c('#4d4d4d',b)

p <- plotAvgProf(tagMatrixList, xlim=c(-1500, 1500),
                 xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency") +
  theme_bw()+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size= 16),
        axis.text = element_text(size= 14),
        legend.position=c(0.1,0.8),legend.title = element_blank(),
        legend.background = element_rect(fill = NA)) +labs(color = 'Group')+
  scale_fill_manual(values = color,aesthetics = "color") ##折线图

Ipaper::write_fig(p,file = "/NAS/yulix/all_G4_data/figure/read_count_frequency.genomic_region.pdf",
                  width = 7,height = 6.5,devices = NULL,res = 300,show = F) 


