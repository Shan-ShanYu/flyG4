rm(list = ls())
Num = "004.1."

#### 原基因tpm ####  
genecounts <- fread("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/genecounts1.txt") %>% as.data.frame
genecounts <- genecounts[,c(1,7:10)]
colnames(genecounts) <- c("Geneid","Female1","Female2","Male1","Male2")
write.table(genecounts,file = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/",Num,"MaleFemale_genecounts.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)

gene_length <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.2.dmel.genelength.txt")

countsTPM <- function(count, efflength){
  RPK <- count/(efflength/1000) ##每千个碱基
  PMSC_rpk <- colSums(RPK)/1e6
  RPK/PMSC_rpk
}
genecounts <- column_to_rownames(genecounts,"Geneid") ##要有行名，不然在下一步计算时会报错
gene_length <- column_to_rownames(gene_length,"gene.name")
gene_tpm <- as.data.frame(apply(genecounts, 2, countsTPM, efflength=gene_length)) ##数据中的每一列应用countsTPM函数，其中efflength为G4len
colnames(gene_tpm) <- colnames(genecounts)

##热图
library(pheatmap)
data2=gene_tpm[which(rowSums(gene_tpm==0)==0),] ##去除为0的值

figure1 <- pheatmap(data2,scale = "row",show_rownames=F,color = colorRampPalette(c("navy", "white", "firebrick3"))(50)) 
##scale = "row"对行进行标准化
##show_rownames=F没有行名

figure1<-last_plot()
ggsave(figure1,filename = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/Picture/",Num,"geneTPM.HeatPlot.pdf"),
       device = "pdf",width = 4,height = 3)  

#### mergeG4 ####
#*mergeG4上下游延伸500bp------------------------------------------------------
rm(list = ls())
merge_g4 <- read.table("/home/yuss/flyG4/result/PQS/001.2.merge.bed",header = T)
merge_g4$center <- round((merge_g4$start+merge_g4$end)/2,0)
merge_g4 <- merge_g4[,-c(2,3)]
merge_g4$start <- merge_g4$center-500
merge_g4$end <- merge_g4$center+500
merge_g4 <- merge_g4[,c(1,13,14,2,3,4,9,10,11)]
write.table(merge_g4,file = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/",Num,"merge.G4wide1K.bed"),
            sep = "\t",row.names = F,col.names = F,quote = F)

rm(list = ls())
Num = "004.1."
#*counts--------------------------------------------------------------------
setwd("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq")
##读取该路径下的所有.txt文件
filename <- list.files("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq",pattern = "^004\\.1\\.SRR166.*\\.txt$")
# for (i in 1:length(filename)){
#   var_name <- gsub('.coverage.txt', '', filename[i]) ##更改文件名
#   assign(var_name, read.table(filename[i], sep = '\t', header = F)) ##assign()函数将一个读取的数据框对象分配给先前定义的变量名 var_name
# }

for (i in 1:length(filename)){
  var_name <- gsub('.coverage.txt', '',filename[i])
  file_data <- read.table(filename[i], sep = '\t', header = F)
  assign(var_name,file_data[,c(2,10)])
}
`004.1.SRR166807_Female_Rep1`$female_rep2 <- `004.1.SRR166808_Female_Rep2`[match(`004.1.SRR166807_Female_Rep1`$V2,`004.1.SRR166808_Female_Rep2`$V2),2]
`004.1.SRR166807_Female_Rep1`$male_rep1 <- `004.1.SRR166809_Male_Rep1`[match(`004.1.SRR166807_Female_Rep1`$V2,`004.1.SRR166809_Male_Rep1`$V2),2]
`004.1.SRR166807_Female_Rep1`$male_rep2 <- `004.1.SRR166810_Male_Rep2`[match(`004.1.SRR166807_Female_Rep1`$V2,`004.1.SRR166810_Male_Rep2`$V2),2]
assign("counts", `004.1.SRR166807_Female_Rep1`[,-1]) #assign() 函数将整理好的表赋值给一个变量counts
colnames(counts)[1] <- "female_rep1"
fwrite(counts,file = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/",Num,"QuantifyG4.counts.txt"),
       sep = "\t",row.names = F,quote = FALSE)

#*差异表达分析--------------------------------------------------------------
rm(list = ls())
Num = "004.1."
##1.构建表达矩阵
counts <- fread("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/004.1.QuantifyG4.counts.txt") %>% as.data.frame()
counts1 <- as.matrix(counts)
pheatmap(log2(counts1+1))

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
library(pheatmap)
#pheatmap(normalized_counts,scale = "row")
pheatmap(log2(normalized_counts+1))

res <- results(dds,contrast = c("group", "male", "female"))  #在差异表达分析中不会对Cook's距离进行截断，即不会自动剔除那些可能影响结果的离群值。
summary(res)
##4.将差异表达结果储存为表格形式
res1 <- data.frame(res,stringsAsFactors = F,check.names = F) #%>% na.omit()#移除含有缺失值（NA）的行
res1$id <- 1:nrow(res1) 
res1 <- na.omit(res1)
# res1 <- merge(res1,G4bed,by.x="id",by.y="id")
# res1[is.na(res1$padj), "padj"]=1 ##is.na(res$padj)p值列中哪些元素是缺失值（NA）,res[is.na(res$padj),"padj"]p值列中是缺失值的那些行，并且在 "padj" 这一列进行操作
# res1[is.na(res1)] = 0 ##所有缺失值（NA）都设置为0
res1$group <- ifelse(res1$padj<0.05&abs(res1$log2FoldChange)>1,ifelse(res1$log2FoldChange>1,"Male-biased","Female-biased"),"Unbias")
table(res1$group)
####和上面代码相同
# res1$group1 <- "Unbias"
# res1[res1$padj < 0.05 & res1$log2FoldChange > 1,"group1"] = "Male-biased"
# res1[res1$padj < 0.05 & res1$log2FoldChange < -1,"group1"] = "Female-biased"
fwrite(res1,file = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/",Num,"QuantifyG4.DEG4.txt"),
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
  annotate(geom="text", x=-7.5, y=70, size=4,label=paste0("Female-biased\n(N=",sum(res1$group=="Female-biased"),")"))+
  annotate(geom="text", x=7, y=70, size=4,label=paste0("Male-biased\n(N=",sum(res1$group=="Male-biased"),")"))+
  scale_color_manual(values = c("#f77aae","#82cbff","grey"))+
  guides(color="none") ##隐藏颜色图例
ggsave(filename = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/Picture/",Num,"DEG4.VolcanoPlot.pdf"),
       device = "pdf",width = 4,height = 3)  

#### kc_all_G4 ####  
#*kc_all_G4上下游延伸500bp----------------------------------------------------
rm(list = ls())
Num = "004.1."
kc_g4 <- read.table("/home/yuss/flyG4/result/PQS/001.2.kc_all.bed",header = T)
kc_g4$center <- round((kc_g4$start+kc_g4$end)/2,0)
kc_g4 <- kc_g4[,-c(2,3)]
kc_g4$start <- kc_g4$center-500
kc_g4$end <- kc_g4$center+500
kc_g4 <- kc_g4[,c(1,13,14,2,3,4,9,10,11)]
write.table(kc_g4,file = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/",Num,"kc_all.G4wide1K.bed"),
            sep = "\t",row.names = F,col.names = F,quote = F)
#*counts---------------------------------------------------------------------
filename <- list.files("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq", pattern = "^004\\.1\\.kcall.*\\.txt$") ##^符号来匹配以"kc_all"开头的文件
for (i in 1:length(filename)) {
  var_name <- gsub('.coverage.txt', '',filename[i])
  file_data <- read.table(filename[i], sep = '\t', header = F)
  assign(var_name,file_data[,c(2,10)])
}
`004.1.kcall_Female_Rep1`$female_rep2 <- `004.1.kcall_Female_Rep2`[match(`004.1.kcall_Female_Rep1`$V2,`004.1.kcall_Female_Rep2`$V2),2]
`004.1.kcall_Female_Rep1`$male_rep1 <- `004.1.kcall_Male_Rep1`[match(`004.1.kcall_Female_Rep1`$V2,`004.1.kcall_Male_Rep1`$V2),2]
`004.1.kcall_Female_Rep1`$male_rep2 <- `004.1.kcall_Male_Rep2`[match(`004.1.kcall_Female_Rep1`$V2,`004.1.kcall_Male_Rep2`$V2),2]
assign("kcall_counts",`004.1.kcall_Female_Rep1`[,-1])
colnames(kcall_counts)[1] <- "female_rep1"
fwrite(kcall_counts,file = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/",Num,"Quantify.kcall_G4.counts.txt"),
       sep = "\t",row.names = F,quote = FALSE)

#*差异表达分析--------------------------------------------------------------------------------------------------
rm(list = ls())
Num = "004.1."
##1.构建表达矩阵
counts <- fread("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/004.1.Quantify.kcall_G4.counts.txt") %>% as.data.frame()
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
res1$group <- ifelse(res1$padj<0.05&abs(res1$log2FoldChange)>1,ifelse(res1$log2FoldChange>1,"Male-biased","Female-biased"),"Unbias")
table(res1$group)
####和上面代码相同
# res1$group1 <- "Unbias"
# res1[res1$padj < 0.05 & res1$log2FoldChange > 1,"group1"] = "Male-biased"
# res1[res1$padj < 0.05 & res1$log2FoldChange < -1,"group1"] = "Female-biased"
fwrite(res1,file = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/",Num,"Quantify.kcall_G4.DEG4.txt"),
       sep = "\t",row.names = F,quote = FALSE)

#### s2_all_G4 ####
rm(list = ls());gc();rm(list = ls())
Num = "004.1."
#*s2_all_G4上下游延伸500bp-------------------------------------------------------------------
s2_g4 <- read.table("/home/yuss/flyG4/result/PQS/001.2.s2_all.bed",header = T)
s2_g4$center <- round((s2_g4$start+s2_g4$end)/2,0)
s2_g4 <- s2_g4[,-c(2,3)]
s2_g4$start <- s2_g4$center-500
s2_g4$end <- s2_g4$center+500
s2_g4 <- s2_g4[,c(1,13,14,2,3,4,9,10,11)]
write.table(s2_g4,file = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/",Num,"s2_all.G4wide1K.bed"),
            sep = "\t",row.names = F,col.names = F,quote = F)
#*counts------------------------------------------------------------------------------------
filename <- list.files("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq", pattern = "^004\\.1\\.s2all.*\\.txt$") ##^符号来匹配以"kc_all"开头的文件
for (i in 1:length(filename)) {
  var_name <- gsub('.coverage.txt', '',filename[i])
  file_data <- read.table(filename[i], sep = '\t', header = F)
  assign(var_name,file_data[,c(2,10)])
}
`004.1.s2all_Female_Rep1`$female_rep2 <- `004.1.s2all_Female_Rep2`[match(`004.1.s2all_Female_Rep1`$V2,`004.1.s2all_Female_Rep2`$V2),2]
`004.1.s2all_Female_Rep1`$male_rep1 <- `004.1.s2all_Male_Rep1`[match(`004.1.s2all_Female_Rep1`$V2,`004.1.s2all_Male_Rep1`$V2),2]
`004.1.s2all_Female_Rep1`$male_rep2 <- `004.1.s2all_Male_Rep2`[match(`004.1.s2all_Female_Rep1`$V2,`004.1.s2all_Male_Rep2`$V2),2]
assign("s2all_counts",`004.1.s2all_Female_Rep1`[,-1])
colnames(s2all_counts)[1] <- "female_rep1"
fwrite(s2all_counts,file = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/",Num,"Quantify.s2all_G4.counts.txt"),
       sep = "\t",row.names = F,quote = FALSE)

#*差异表达分析--------------------------------------------------------------------------------------------------
##1.构建表达矩阵
Num = "004.1."
counts <- fread("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/004.1.Quantify.s2all_G4.counts.txt") %>% as.data.frame()
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
res1$group <- ifelse(res1$padj<0.05&abs(res1$log2FoldChange)>1,ifelse(res1$log2FoldChange>1,"Male-biased","Female-biased"),"Unbias")
table(res1$group)
# Female-biased   Male-biased        Unbias 
# 1594           994          6790 
####和上面代码相同
# res1$group1 <- "Unbias"
# res1[res1$padj < 0.05 & res1$log2FoldChange > 1,"group1"] = "Male-biased"
# res1[res1$padj < 0.05 & res1$log2FoldChange < -1,"group1"] = "Female-biased"
fwrite(res1,file = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/",Num,"Quantify.s2all_G4.DEG4.txt"),
       sep = "\t",row.names = F,quote = FALSE)


#### overlapG4 ####
#*overlapG4上下游扩500bp---------------------------------------------------------
rm(list = ls())
Num = "004.1."
overlapg4 <- read.table("/home/yuss/flyG4/result/PQS/001.2.overlap.bed",header = T)
overlapg4$center <- round((overlapg4$start+overlapg4$end)/2,0)
overlapg4 <- overlapg4[,-c(2,3)]
overlapg4$start <- overlapg4$center-500
overlapg4$end <- overlapg4$center+500
overlapg4 <- overlapg4[,c(1,13,14,2,3,4,9,10,11)]
write.table(overlapg4,file = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/",Num,"overlap.G4wide1K.bed"),
            sep = "\t",row.names = F,col.names = F,quote = F)

#*counts------------------------------------------------------------------------------------
filename <- list.files("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq", pattern = "^004\\.1\\.overlap.*\\.txt$") ##^符号来匹配以"kc_all"开头的文件
for (i in 1:length(filename)) {
  var_name <- gsub('.coverage.txt', '',filename[i])
  file_data <- read.table(filename[i], sep = '\t', header = F)
  assign(var_name,file_data[,c(2,10)])
}
`004.1.overlap_Female_Rep1`$female_rep2 <- `004.1.overlap_Female_Rep2`[match(`004.1.overlap_Female_Rep1`$V2,`004.1.overlap_Female_Rep2`$V2),2]
`004.1.overlap_Female_Rep1`$male_rep1 <- `004.1.overlap_Male_Rep1`[match(`004.1.overlap_Female_Rep1`$V2,`004.1.overlap_Male_Rep1`$V2),2]
`004.1.overlap_Female_Rep1`$male_rep2 <- `004.1.overlap_Male_Rep2`[match(`004.1.overlap_Female_Rep1`$V2,`004.1.overlap_Male_Rep2`$V2),2]
assign("overlap_counts",`004.1.overlap_Female_Rep1`[,-1])
colnames(overlap_counts)[1] <- "female_rep1"
fwrite(overlap_counts,file = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/",Num,"Quantify.overlap_G4.counts.txt"),
       sep = "\t",row.names = F,quote = FALSE)

#*差异表达分析--------------------------------------------------------------------------------------------------
##1.构建表达矩阵
counts <- fread("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/004.1.Quantify.overlap_G4.counts.txt") %>% as.data.frame()
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
res1$group <- ifelse(res1$padj<0.05&abs(res1$log2FoldChange)>1,ifelse(res1$log2FoldChange>1,"Male-biased","Female-biased"),"Unbias")
table(res1$group)
# Female-biased   Male-biased        Unbias 
# 1304           819          5548 
# res1$group <- ifelse(res1$pvalue<0.05&abs(res1$log2FoldChange)>1,ifelse(res1$log2FoldChange>1,"Male-biased","Female-biased"),"Unbias")
# table(res1$group)
# Female-biased   Male-biased        Unbias 
# 1558          1199          4914
# res1$group <- ifelse(res1$pvalue<0.05&abs(res1$log2FoldChange)>2,ifelse(res1$log2FoldChange>2,"Male-biased","Female-biased"),"Unbias")
# table(res1$group)
# Female-biased   Male-biased        Unbias 
# 789           701          6181 
# res1$group <- ifelse(res1$pvalue<0.05&abs(res1$log2FoldChange)>0.5,ifelse(res1$log2FoldChange>0.5,"Male-biased","Female-biased"),"Unbias")
# table(res1$group)
# Female-biased   Male-biased        Unbias 
# 1634          1239          4798
####和上面代码相同
# res1$group1 <- "Unbias"
# res1[res1$padj < 0.05 & res1$log2FoldChange > 1,"group1"] = "Male-biased"
# res1[res1$padj < 0.05 & res1$log2FoldChange < -1,"group1"] = "Female-biased"
fwrite(res1,file = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/",Num,"Quantify.overlap_G4.DEG4.txt"),
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
ggsave(filename = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/Picture/",Num,"OverlapDEG4.VolcanoPlot.pdf"),
       device = "pdf",width = 4,height = 3)  
