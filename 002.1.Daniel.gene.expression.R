rm(list = ls());gc();rm(list = ls())#清空
Num = "002.1."
##BiocManager::install("edgeR")
library(edgeR)
#### 1.差异表达分析 ####
#导入数据
rawdata <- fread('/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/counts.matrix2') %>% as.data.frame()
colnames(rawdata) <- c('Geneid','kc','s2')
counts <- column_to_rownames(rawdata,var = "Geneid")
write.table(counts,file = paste0("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/",Num,"counts.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)
group <- rep(c('kc','s2'),each=1)
#1.构建DGEList对象
dgelist <- DGEList(counts = counts, group = group) #参数很简单 你的矩阵和分组文件
#2.数据过滤
#原表达矩阵基因数太大，可能存在某些基因没有表达，预先过滤
#保留在至少在一个样本里有表达的基因(CPM > 1)
keep <- rowSums(cpm(dgelist)>1) >= 1
dgelist <- dgelist[keep,,keep.lib.sizes=FALSE]
#3.标准化
#考虑到测序深度不同, 我们需要对其进行标准化, 避免文库大小不同导致的分析误差
dgelist <- calcNormFactors(dgelist, method = 'TMM')
dgelist$samples
#4.差异表达分析
#不同差异表达分析工具的目标就是预测出dispersion(离散值), 有了离散值就能够计算p值
dgelist.bcv <- dgelist
bcv <- 0.4
et <- exactTest(dgelist.bcv, dispersion = bcv ^ 2)
#使用decideTestsDGE看下有多少基因上调, 多少基因下调
#设置p.value=0.05
gene1 <- decideTestsDGE(et, p.value = 0.05, lfc = 0)
summary(gene1)
head(gene1)
colnames(gene1) <- 'type'
library(ggplot2)
#差异基因结果
results <- cbind(dgelist$counts,et$table,gene1)
results <- rownames_to_column(results, 'gene_id')
results$siginifi <- ifelse(abs(results$type) == 1,
                           ifelse(results$type == 1, "up", "down"),"nosig") #UP：差异显著且上调的基因。DOWN：差异显著且下调的基因。NOT：差异不显著的基因。
table(results$siginifi)
colnames(results)
results_select <- subset(results, abs(type)==1)
write.table(results,file = paste0("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/",Num,"DE.s2vskc.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)
#5.火山图可视化
library(ggplot2)
ggplot(results,aes(x=logFC,y=-log10(PValue))) +#黑白主题，base_size可以设置某种主题内基本字体的大小
  geom_point(aes(color=siginifi),alpha=0.5,size=2) +
  theme_bw(base_size = 16)+
  theme(aspect.ratio = 1,
        #plot.title = element_text(hjust = 0.5), # #设置标题居中
        axis.text = element_text(color = 'black'),
        axis.title = element_text(color = 'black')) +
  scale_color_manual(name = '',
                     values = c('up'='#D6604D','nosig'='grey','down'='#74ADD1'), #手动设置颜色时调整颜色的因子顺序
                     label = c('up'='S2','noSig'='nosig','down'='KC'))+
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',size = 0.8) +
  geom_vline(xintercept = c(-1,1),lty = 'dashed',size = 0.8) +
  #scale_x_continuous(breaks = c(-6,-4,-2,-1,0,1,2,4,6)) +
  ## 修改坐标轴
  labs(x="log2 Fold Change",y="-Log10 (p-value)")

#ggsave(file = "up_GI3.png",p1,width = 10,height = 8)

#### 2.gene与kcG4 s2G4交集后合并表达量 ####
#读取交集的表(pqs和s2、kc peak取交集得到哪些pqs可以形成G4，形成G4的bed（s2_all.bed）再与gene.bed取交集得到gene上有S2的G4)
gene.kc <- fread('/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.gene.kc.bed') %>% as.data.frame()
gene.s2 <- fread('/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.gene.s2.bed') %>% as.data.frame()
results <- fread('/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.DE.s2vskc.txt')
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
write.table(gene.kc.s2,file = paste0("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/",Num,"gene.kc.s2.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)

gene.kcG4 <- gene.kc.s2[gene.kc.s2$sum==1,]
gene.s2G4 <- gene.kc.s2[gene.kc.s2$sum==2,]
gene.kcs2G4 <- gene.kc.s2[gene.kc.s2$sum==3,]
gene.nonG4 <- gene.kc.s2[gene.kc.s2$sum==0,]
#按行合并
gene.kcG4$type <- 'gene.kcG4'
gene.s2G4$type <- 'gene.s2G4'
gene.kcs2G4$type <- 'gene.kcs2G4'
gene.nonG4$type <- 'gene.nonG4'
gene.merge.all <- bind_rows(gene.kcs2G4,gene.kcG4,gene.s2G4,gene.nonG4)
#合并log2Foldchange
gene.merge.all$logFC <- results[match(gene.merge.all$id,results$gene_id),4]
#RNAseq数据过滤这一步过滤了一万左右表达量低的reads,导致gene.merge.all与resault合并时，gene.merge.all中的logFC匹配不到，有NA值，这要去除掉
gene.merge.all=subset(gene.merge.all, logFC!= 'NA') #只去除logFC列含有NA值的行
gene.merge.all$type <- factor(gene.merge.all$type,levels = c("gene.nonG4","gene.kcG4","gene.s2G4","gene.kcs2G4"))
#可视化含有G4的基因和不含有G4的基因的表达水平
my_comparisons = list(c("gene.nonG4","gene.kcG4"),c("gene.kcG4","gene.s2G4"),
                      c("gene.s2G4","gene.kcs2G4"))
b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(5)
color <- c('#4d4d4d',b)
library(ggpubr) #是一个基于ggplot2的计算工具包
ggplot(data = gene.merge.all,aes(x=type,y=logFC,color=type)) +
  geom_boxplot() +
  geom_point()+
  geom_jitter(width = 0.25) + #点是抖动的且没有超出箱线图的宽度
  scale_color_manual(values = color)+ #更改箱线图箱子的颜色
  theme_bw() +
  stat_compare_means(comparisons = my_comparisons,label.y = c(rep(c(16.5,18),3)),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4) +
  ylab("log2 Fold Change") +
  labs(title="S2 vs KC") +
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") + #没有图例
  coord_cartesian(ylim = c(-16,19))

#### 3.上下调表达的基因中含有G4的数量，看一下G4是否会影响基因的表达 ####
#把results$siginifi合并到gene.merge.all中以便统计四种类型G4的数量，之前只合并logFC
gene.merge.all$siginifi <- results[match(gene.merge.all$id,results$gene_id),8]
table(gene.merge.all$type,gene.merge.all$siginifi)
gene.merge.all.table <- table(gene.merge.all$type,gene.merge.all$siginifi) %>% as.data.frame()
#使用spread函数将长数据转换为宽数据
gene.merge.all.table1 <- spread(gene.merge.all.table,Var2,Freq)
colnames(gene.merge.all.table1)[1] <- "type"
b = colorRampPalette(colors = c('#fddbc7', '#67001f'))(5)
color <- c('#4d4d4d',b)
#在s2细胞中上调表达的基因是否含有G4的数量
library(ggplot2)
ggplot(gene.merge.all.table1, aes(x=type,y=up,fill=type)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(0.7),width = 0.7,color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值
  coord_cartesian(ylim = c(0,500)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = color) + ##更改颜色
  geom_text(aes(label=up,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题(左下边框，没有网格线)
  ylab("Number of S2 G4-containing genes") + ##通过bquote函数给图标签添加上下标
  theme(axis.title.y = element_text(size = 16), ##y坐标轴标题字体大小
        axis.title.x = element_blank(), ##删除x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.position = "none") 
#在kc细胞中上调表达的基因是否含有G4的数量
ggplot(gene.merge.all.table1, aes(x=type,y=down,fill=type)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(0.7),width = 0.7,color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值
  coord_cartesian(ylim = c(0,400)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = color) + ##更改颜色
  geom_text(aes(label=down,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题(左下边框，没有网格线)
  ylab("Number of KC G4-containing genes") + ##通过bquote函数给图标签添加上下标
  theme(axis.title.y = element_text(size = 16), ##y坐标轴标题字体大小
        axis.title.x = element_blank(), ##删除x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.position = "none") 

#### 4.promoter与kcG4 s2G4交集后合并表达量####
#提取promoter.bed
gene.sort <- fread('/home/yuss/flyG4/data/ref/gene.sort.bed')
gene.sort.p <- gene.sort[gene.sort$V5=='+',] ##提取正链的信息
gene.sort.m <- gene.sort[gene.sort$V5=='-',]
promoter.p <- data.frame(gene.sort.p$V1,gene.sort.p$V2-1500,gene.sort.p$V2+1500,gene.sort.p$V4,gene.sort.p$V5)
colnames(promoter.p) <- c("chr","start","end","gene_id","strand")
promoter.m <- data.frame(gene.sort.m$V1,gene.sort.m$V3-1500,gene.sort.m$V3+1500,gene.sort.m$V4,gene.sort.m$V5)
colnames(promoter.m) <- c("chr","start","end","gene_id","strand")
promoter <- rbind(promoter.p,promoter.m)
promoter$start <- ifelse(promoter$start <0, 0, promoter$start)
write.table(promoter,file = '/home/yuss/flyG4/data/ref/promoter.bed',
            sep = '\t',col.names = T,row.names = F,quote = F)
#读取promoter与s2,kc G4交集的文件
promoter.kc <- fread('/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.promoter.kc.bed')
promoter.s2 <- fread('/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.promoter.s2.bed')
#合并两个表
promoter.kc.s2 <- bind_cols(promoter.kc,promoter.s2$V6)
#修改列名,合并表达量
colnames(promoter.kc.s2) <- c("chr","start","end","gene_id","strand","kc","s2")
results <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.DE.s2vskc.txt") %>% data.frame()
promoter.kc.s2$kc <- ifelse(promoter.kc.s2$kc==0,0,1)
promoter.kc.s2$s2 <- ifelse(promoter.kc.s2$s2==0,0,1)
promoter.kc.s2$logFC <- results[match(promoter.kc.s2$gene_id,results$gene_id),4]
promoter.kc.s2 <- subset(promoter.kc.s2, logFC!= 'NA')
promoter.kc.s2$logFC <- as.numeric(promoter.kc.s2$logFC)
promoter.kc.s2$logFC.round <- round(promoter.kc.s2$logFC, 0)
#计数
promoter.kc.logFC.num <- data.frame()
promoter.kc.logFC.num <- table(promoter.kc.s2$logFC.round, promoter.kc.s2$kc) %>% as.data.frame()
promoter.s2.logFC.num <- data.frame()
promoter.s2.logFC.num <- table(promoter.kc.s2$logFC.round,promoter.kc.s2$s2) %>% as.data.frame()
#使用spread函数将gd1_long长数据转换为宽数据gd1_wide
promoter.kc.logFC.num <- spread(promoter.kc.logFC.num,Var2,Freq) ##长数据转换为宽数据，Var2为需要分解的变量，Freq为分解后的列的取值
colnames(promoter.kc.logFC.num) <- c("logFC","non","have")
promoter.kc.logFC.num$non %<>% as.numeric(.)
promoter.kc.logFC.num$have %<>% as.numeric(.)
promoter.s2.logFC.num <- spread(promoter.s2.logFC.num,Var2,Freq)
colnames(promoter.s2.logFC.num) <- c("logFC","non","have")
promoter.s2.logFC.num$non %<>% as.numeric(.)
promoter.s2.logFC.num$have %<>% as.numeric(.)
#计算不同log2FC上promoter含有G4的比例（有G4/有G4+无G4）
promoter.kc.logFC.num$sum <- rowSums(promoter.kc.logFC.num[,c(2,3)])
promoter.kc.logFC.num$ratio <- promoter.kc.logFC.num$have/promoter.kc.logFC.num$sum 
promoter.kc.logFC.num$ratio <- round(promoter.kc.logFC.num$ratio, 2) ##保留两位小数
promoter.s2.logFC.num$ratio <- promoter.s2.logFC.num$have/rowSums(promoter.s2.logFC.num[,c(2,3)])
promoter.s2.logFC.num$ratio <- round(promoter.s2.logFC.num$ratio,2)
#删除log2FC=0这一行
promoter.kc.logFC.num <- promoter.kc.logFC.num[-14,]
promoter.s2.logFC.num <- promoter.s2.logFC.num[-14,]
#可视化含有G4的promoter和不含有G4的promoter的表达水平
#增加type这一列后面柱形图分类
promoter.kc.logFC.num$num <- c(-13:-1,1:16)
promoter.kc.logFC.num$type <- ifelse(promoter.kc.logFC.num$num>0,'s2','kc')
promoter.s2.logFC.num$num <- c(-13:-1,1:16)
promoter.s2.logFC.num$type <- ifelse(promoter.s2.logFC.num$num>0,'s2','kc')
ggplot(promoter.kc.logFC.num, aes(x=logFC,y=ratio*100,fill=type)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(),width = 0.9,
           color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值，#width设置矩形条的宽度
  scale_fill_manual(values = c('s2'='#D6604D','kc'='#74ADD1')) +#手动设置颜色时调整颜色的因子顺序
  coord_cartesian(ylim = c(0,50)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) + ##消除x轴与绘图区的间隙
  #geom_text(aes(label=Freq,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题
  ylab("% G4 genes in KC cells") + 
  #ylab("Normalized density (No./kb/GC%)") +
  xlab("log2 Fold Change") +
  #scale_fill_discrete(name="group") + ##修改图例标题名称
  theme(axis.title.y = element_text(size = 16), ##y坐标轴标题字体大小
        axis.title.x = element_text(size = 16), ##x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.position = "none") ##删除图例

ggplot(promoter.s2.logFC.num, aes(x=logFC,y=ratio*100,fill=type)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(),width = 0.9,
           color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值，#width设置矩形条的宽度
  scale_fill_manual(values = c('s2'='#D6604D','kc'='#74ADD1')) +#手动设置颜色时调整颜色的因子顺序
  coord_cartesian(ylim = c(0,30)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) + ##消除x轴与绘图区的间隙
  #geom_text(aes(label=Freq,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题
  ylab("% G4 genes in S2 cells") + 
  #ylab("Normalized density (No./kb/GC%)") +
  xlab("log2 Fold Change") +
  #scale_fill_discrete(name="group") + ##修改图例标题名称
  theme(axis.title.y = element_text(size = 16), ##y坐标轴标题字体大小
        axis.title.x = element_text(size = 16), ##x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.position = "none") ##删除图例标题名称

ggplot(promoter.kc.logFC.num, aes(x=logFC,y=have,fill=type)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(),width = 0.9,
           color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值，#width设置矩形条的宽度
  scale_fill_manual(values = c('s2'='#D6604D','kc'='#74ADD1')) +#手动设置颜色时调整颜色的因子顺序
  coord_cartesian(ylim = c(0,300)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) + ##消除x轴与绘图区的间隙
  #geom_text(aes(label=Freq,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题
  ylab("Number of G4 genes in KC cells") + 
  #ylab("Normalized density (No./kb/GC%)") +
  xlab("log2 Fold Change") +
  #scale_fill_discrete(name="group") + ##修改图例标题名称
  theme(axis.title.y = element_text(size = 16), ##y坐标轴标题字体大小
        axis.title.x = element_text(size = 16), ##x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小) 
        legend.title = element_blank(),
        legend.position = c(.90, .90),
        legend.text=element_text(size=14)) 

ggplot(promoter.s2.logFC.num, aes(x=logFC,y=have,fill=type)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(),width = 0.9,
           color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值，#width设置矩形条的宽度
  scale_fill_manual(values = c('s2'='#D6604D','kc'='#74ADD1')) +#手动设置颜色时调整颜色的因子顺序
  coord_cartesian(ylim = c(0,300)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) + ##消除x轴与绘图区的间隙
  #geom_text(aes(label=Freq,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题
  ylab("Number of G4 genes in S2 cells") + 
  #ylab("Normalized density (No./kb/GC%)") +
  xlab("log2 Fold Change") +
  #scale_fill_discrete(name="group") + ##修改图例标题名称
  theme(axis.title.y = element_text(size = 16), ##y坐标轴标题字体大小
        axis.title.x = element_text(size = 16), ##x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.title = element_blank(), ##删除图列标题
        legend.position = c(.90, .90),
        legend.text=element_text(size=14))  ##图例字体大小

#### 5.GO 富集分析 ####
rm(list = ls());gc()
library(ChIPseeker)
library(clusterProfiler)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
path <- "/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed"
files <- dir(path)
filepath <- sapply(files, function(x){
  paste(path,x,sep="/")
})
peak=list(non_eG4=filepath[[3]],kc_specific=filepath[[1]],s2_specific=filepath[[5]],
          overlap=filepath[[4]],merge=filepath[[2]])
peakAnnoList = list()
for (i in 1:length(peak)) {
  peakAnnoList[[i]] <- lapply(peak[[i]], annotatePeak, TxDb=txdb,
                             tssRegion=c(-500,500),verbose=FALSE)
}
gene=list()
for (i in 1:length(peak)) {
  gene[[i]] <- lapply(peakAnnoList[[i]], function(i) as.data.frame(i)$geneId)
}
list=list()
for (i in 1:length(peak)) {
  list[[i]] <- data.frame(t(data.frame(gene[[i]])))
  list[[i]] <- unlist(as.character(list[[i]])) #unlist() 函数用于将列表转换为向量
}
names(list) <- c("non_eG4","kc_specific","s2_specific","overlap","merge")

#id转化
library(org.Dm.eg.db)
columns(org.Dm.eg.db)
gene.bitr <- list()
for (i in c(1:5)) {
  gene.bitr[[i]] = bitr(geneID = list[[i]],
                          fromType = "FLYBASE",
                          toType = c("ENTREZID","SYMBOL"),
                          OrgDb = "org.Dm.eg.db")
}
gene.ENTREZID <- list() #取$ENTREZID
for (i in c(1:5)) {
  gene.ENTREZID[[i]] = gene.bitr[[i]]$ENTREZID
}
#GO富集
enrich.go.BP <- list()
for(i in c(1:5)) {
  enrich.go.BP[[i]] =enrichGO(gene = gene.ENTREZID[[i]],
                              OrgDb = org.Dm.eg.db,
                              keyType = "ENTREZID",
                              ont = "BP",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.1)
}
#批量导出GO BP结果
outpath <-'/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/'
out_filename <- sapply(names(list), function(x){
  paste(x,'.csv',sep='')
})
outfilepath <- sapply(out_filename,function(x){
  paste0(outpath,Num,x)
})
for (i in 1:length(enrich.go.BP)) {
  write.table(enrich.go.BP[[i]],file = outfilepath[[i]],
              sep = '\t',row.names = FALSE,quote = FALSE)
}

a <- names(list)
data = enrich.go.BP
for (i in 1:length(data)) {
  data[[i]] <- arrange(data[[i]],p.adjust)
  data[[i]] <- data[[i]][1:7,]
  data[[i]] <- data[[i]][,c("Description","p.adjust")]
  data[[i]]$group <- a[i]
}
df <- do.call('rbind',data)
df <- dplyr::mutate(df,p.value=-log10(p.adjust)) #新增一列
df <- df[,c(1,3,4)]
max(df$p.value)
min(df$p.value)
df <- spread(data = df, key = group, value = p.value) #长数据转宽数据
df[is.na(df)] <- 0
library(Hmisc)
df$Description <- capitalize(df$Description) #只是将首字母修改为大写
df <- as.data.frame(column_to_rownames(df,var="Description"))
#G4的GO富集（merge）
# data_merge <- enrich.go.BP[[5]]
# data_merge <- arrange(data_merge, p.adjust)
# data_merge <- data_merge[,c("Description","p.adjust")]
dotplot(enrich.go.BP[[5]], showCategory=10,title="EnrichmentGO_Merge_BP") #与细胞形态发生有关
data_kc_specific <- enrich.go.BP[[2]]
data_kc_specific <- arrange(data_kc_specific, p.adjust)
data_kc_specific <- data_kc_specific[,c("Description","p.adjust")]
data_non_eG4 <- arrange(enrich.go.BP[[1]],p.adjust)[,c("Description","p.adjust")]
data_s2_specific <- arrange(enrich.go.BP[[3]],p.adjust)[,c("Description","p.adjust")]
dotplot(enrich.go.BP[[1]], showCategory=10,title="EnrichmentGO_non_eG4_BP",orderBy="p.adjust")
dotplot(enrich.go.BP[[2]], showCategory=10,title="EnrichmentGO_kc_specific_BP",orderBy="p.adjust") 
dotplot(enrich.go.BP[[3]], showCategory=10,title="EnrichmentGO_s2_specific_BP",orderBy="p.adjust")
dotplot(enrich.go.BP[[4]], showCategory=10,title="EnrichmentGO_overlap_BP",orderBy="p.adjust")

# BiocManager::install("rGREAT")
library(rGREAT)
# library(ChIPseeker)
# set.seed(123)
# kc.all.peak <- readPeakFile('/home/yuss/flyG4/result/PQS/001.2.kc_all.bed')
# a <- fread('/home/yuss/flyG4/result/PQS/001.2.kc_all.bed')
# res <- great(a, "GO:BP", "TxDb.Dmelanogaster.UCSC.dm6.ensGene")



gene.sort <- fread('/home/yuss/flyG4/data/ref/gene.sort.bed')
gene.sort.p <- gene.sort[gene.sort$V5=='+',] ##提取正链的信息
gene.sort.m <- gene.sort[gene.sort$V5=='-',]
promoter.p <- data.frame(gene.sort.p$V1,gene.sort.p$V2-100,gene.sort.p$V2+100,gene.sort.p$V4,gene.sort.p$V5,gene.sort.p$V6)
colnames(promoter.p) <- c("chr","start","end","gene_id","strand","name")
promoter.m <- data.frame(gene.sort.m$V1,gene.sort.m$V3-100,gene.sort.m$V3+100,gene.sort.m$V4,gene.sort.m$V5,gene.sort.m$V6)
colnames(promoter.m) <- c("chr","start","end","gene_id","strand","name")
promoter <- rbind(promoter.p,promoter.m)
promoter$start <- ifelse(promoter$start <0, 0, promoter$start)
write.table(promoter,file = '/home/yuss/flyG4/data/ref/TSS10.bed',
            sep = '\t',col.names = T,row.names = F,quote = F)
