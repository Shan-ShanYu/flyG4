rm(list = ls())
Num = "007.5."
kcg4 <- read.table("/home/yuss/flyG4/result/PQS/001.2.kc_specific.bed",header = T)
kcg4$center <- round((kcg4$start+kcg4$end)/2,0)
kcg4 <- kcg4[,-c(2,3)]
kcg4$start <- kcg4$center-500
kcg4$end <- kcg4$center+500
kcg4 <- kcg4[,c(1,13,14,2,3,4,9,10,11)]
write.table(kcg4,file = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/",Num,"Kcspecific.G4wide1K.bed"),
            sep = "\t",row.names = F,col.names = F,quote = F)

s2g4 <- read.table("/home/yuss/flyG4/result/PQS/001.2.s2_specific.bed",header = T)
s2g4$center <- round((s2g4$start+s2g4$end)/2,0)
s2g4 <- s2g4[,-c(2,3)]
s2g4$start <- s2g4$center-500
s2g4$end <- s2g4$center+500
s2g4 <- s2g4[,c(1,13,14,2,3,4,9,10,11)]
write.table(s2g4,file = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/",Num,"S2specific.G4wide1K.bed"),
            sep = "\t",row.names = F,col.names = F,quote = F)
#### overlap kcs2RNAseq####
#*counts------------------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())
Num = "007.5."
getwd()
setwd("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/coverage")
filename <- list.files("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/coverage", pattern = "^007\\.5\\.overlap.*\\.txt$")
for (i in 1:length(filename)){
  var_name <- gsub('\\.coverage\\.txt$|RNAseq_|007\\.5\\.', '', filename[i])
  file_data <- read.table(filename[i], sep = '\t', header = F)
  assign(var_name,file_data[,c(2,10)])
}

`overlap_kcRep1`$kcrep2 <- `overlap_kcRep2`[match(`overlap_kcRep1`$V2,`overlap_kcRep2`$V2),2]
`overlap_kcRep1`$s2rep1 <- `overlap_s2Rep1`[match(`overlap_kcRep1`$V2,`overlap_s2Rep1`$V2),2]
`overlap_kcRep1`$s2rep2 <- `overlap_s2Rep2`[match(`overlap_kcRep1`$V2,`overlap_s2Rep2`$V2),2]
assign("overlap_counts",`overlap_kcRep1`[,-1])
colnames(overlap_counts)[1] <- "kcrep1"
fwrite(overlap_counts,file = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/",Num,"Quantify.overlap_G4.counts.txt"),
       sep = "\t",row.names = F,quote = FALSE)

#*差异表达分析--------------------------------------------------------------------------------------------------
##1.构建表达矩阵
counts <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.5.Quantify.overlap_G4.counts.txt") %>% as.data.frame()
coldata <- data.frame(row.names = colnames(counts),
                      group = c("female","female","male","male"))
coldata$group <- as.factor(coldata$group)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ group)
##2.dds标准化
dds = DESeq(dds)
##3.获取标准化后的数据
normalized_counts <- counts(dds,normalized=T)
res <- results(dds,contrast = c("group", "male", "female"))  #在差异表达分析中不会对Cook's距离进行截断，即不会自动剔除那些可能影响结果的离群值。
summary(res)
##4.将差异表达结果储存为表格形式
res1 <- data.frame(res,stringsAsFactors = F,check.names = F) #%>% na.omit()#移除含有缺失值（NA）的行
res1$id <- 1:nrow(res1) 
res1 <- na.omit(res1)
# res1 <- merge(res1,G4bed,by.x="id",by.y="id")
# res1[is.na(res1$padj), "padj"]=1 ##is.na(res$padj)p值列中哪些元素是缺失值（NA）,res[is.na(res$padj),"padj"]p值列中是缺失值的那些行，并且在 "padj" 这一列进行操作
# res1[is.na(res1)] = 0 ##所有缺失值（NA）都设置为0
res1$group <- ifelse(res1$pvalue<0.05&abs(res1$log2FoldChange)>1,ifelse(res1$log2FoldChange>1,"Male-biased","Female-biased"),"Unbias")
table(res1$group)
# Female-biased   Male-biased        Unbias 
# 508           555          2291 
####和上面代码相同
# res1$group1 <- "Unbias"
# res1[res1$padj < 0.05 & res1$log2FoldChange > 1,"group1"] = "Male-biased"
# res1[res1$padj < 0.05 & res1$log2FoldChange < -1,"group1"] = "Female-biased"
fwrite(res1,file = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/",Num,"Quantify.overlap_G4.DEG4.txt"),
       sep = "\t",row.names = F,quote = FALSE)
##5.可视化
library(ggplot2)
ggplot(res1,aes(log2FoldChange,-log10(padj),color=group))+
  geom_point()+
  theme_bw()+
  xlab(bquote(Log[2]~Fold~Change))+
  ylab(bquote(Log[10]~P~Value))+
  theme(
    axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12))+
  geom_vline(xintercept = -1,linetype="longdash")+
  geom_vline(xintercept = 1,linetype="longdash")+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  annotate(geom="text", x=-8.7, y=75, size=4,label=paste0("Female-biased\n(N=",sum(res1$group=="Female-biased"),")"))+
  annotate(geom="text", x=7.5, y=75, size=4,label=paste0("Male-biased\n(N=",sum(res1$group=="Male-biased"),")"))+
  scale_color_manual(values = c("#f77aae","#82cbff","grey"))+
  coord_cartesian(ylim = c(0,80))+
  guides(color="none") ##隐藏颜色图例
# ggsave(filename = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/Picture/",Num,"OverlapDEG4.VolcanoPlot.pdf"),
       device = "pdf",width = 4,height = 3)  



DEG4 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.5.Quantify.overlap_G4.DEG4.txt") %>% as.data.frame()
G4bed <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.5.overlap.G4wide1K.bed") %>% as.data.frame()
G4bed$id <- 1:nrow(G4bed)
DEG4 <- merge(DEG4,G4bed,by.x="id",by.y="id")
DEG4 %<>% dplyr::filter(!V1 %in% c(4,"Y")) ##V1列的元素不在向量 (4, "Y") 中的行将被保留
DEG4$group %<>% gsub("-biased","",.) ##其中的 . 表示当前正在处理的元素

#*Distribution分布-----------------------------------------------------------------------
df = table(DEG4$group,DEG4$V1) %>% as.matrix() 
df
df <- apply(df, 1, function(x)x/sum(x)*100) %>% as.data.frame()
df$chr = row.names(df)
df %<>% pivot_longer(cols = 1:3) #将df的列1、2、3进行转换，使其从宽格式变为长格式
ggplot(df, aes(chr,value,fill=name))+
  geom_bar(stat="identity",position = "dodge") + #position = "dodge"：这表示将柱状图进行分组显示，使得不同组的柱子并排排列,即绘制分组柱状图
  theme_bw() +
  scale_fill_manual(values = c("#f77aae","#82cbff","grey")) +
  scale_y_continuous(expand=c(0,0),limits = c(0,35))+ #scale_y_continuous() 可以用于修改数据本身的刻度范围，而 coord_cartesian() 则是调整图形的可视化范围，不会改变数据。
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12),
        legend.position = c(0.5, 0.9),
        legend.direction = "horizontal",
        legend.title = element_blank())+
  ylab("Proportion of G4s (%)") +
  xlab("")+
  labs(title="Overlap eG4") 
ggsave(filename = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/Picture/",Num,"Overlap.DEG4.DistributionPlot.pdf"),
       device = "pdf",width = 4,height = 3)