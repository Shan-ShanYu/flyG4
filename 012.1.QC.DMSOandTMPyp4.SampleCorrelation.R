#*质检，查看样本是否离群-----------------------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())#清空
Num = "012.1."
#### 1.counts ####
summary_readcounts <- fread("/home/yuss/flyG4/result/qians.RNAseq/merged_counts.txt") %>% as.data.frame()
readcounts <- summary_readcounts[7:ncol(summary_readcounts)]
colnames(readcounts) <- unlist(lapply(strsplit(colnames(readcounts), "[/]"),"[[",8))
# strsplit(colnames(readcounts), "[/]") 使用 strsplit 函数，根据斜杠字符 ("/") 将每个列名分割成一个字符串向量。这将生成一个列表，列表中的每个元素是一个字符串向量，包含分割后的各部分。
# lapply(..., "[[", 8) 使用 lapply 函数来遍历上一步生成的列表，并提取每个字符串向量的第 8 个元素，这是分割后的文件名部分。
# unlist(...) 将前一步得到的提取文件名的结果合并成一个字符向量。
colnames(readcounts) <- unlist(lapply(strsplit(colnames(readcounts),"[.]"),"[[",1))
colnames(readcounts) <- sub("mix_", "", colnames(readcounts))
colnames(readcounts) <- sub("_trim1", "", colnames(readcounts))

# 102_50_3	TMP50rep3
# 114_25_3	TMP25rep3
# 126_C1	DMSOrep1
# 150_100_1	TMP100rep1
# 162_50_2	TMP50rep2
# 174_100_3	TMP100rep3
# 18_C2	DMSOrep2
# 30_100_2	TMP100rep2
# 54_C3	DMSOrep3
# 6_50_1	TMP50rep1
# 66_25_2	TMP25rep2
# 78_25_1	TMP25rep1
sample_mapping <- c("102_50_3" = "TMP_50_rep3",
                    "114_25_3" = "TMP_25_rep3",
                    "126_C1"   = "DMSO_rep1",
                    "150_100_1" = "TMP_100_rep1",
                    "162_50_2"  = "TMP_50_rep2",
                    "174_100_3" = "TMP_100_rep3",
                    "18_C2"     = "DMSO_rep2",
                    "30_100_2"  = "TMP_100_rep2",
                    "54_C3"     = "DMSO_rep3",
                    "6_50_1"    = "TMP_50_rep1",
                    "66_25_2"   = "TMP_25_rep2",
                    "78_25_1"   = "TMP_25_rep1")
# 提取旧列名
old_colnames <- colnames(readcounts)

# 将旧列名映射到新列名
new_colnames <- sample_mapping[old_colnames]

# 使用新列名替换旧列名
colnames(readcounts) <- new_colnames
row.names(readcounts) <- summary_readcounts$Geneid
readcounts<-readcounts[apply(readcounts, 1, sum)>0,]

sorted_colnames <- colnames(readcounts)[order(colnames(readcounts))]
sorted_counts <- readcounts[, sorted_colnames]

sorted_counts$geneid <- rownames(sorted_counts)
write.table(sorted_counts,file = paste0("/home/yuss/flyG4/result/qians.RNAseq/",Num,"counts.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)

rm(list = ls());gc();rm(list = ls())#清空
Num = "012.1."
#### 2.计算TPM ####
#之前计算过基因长度
gene.length <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.2.dmel.genelength.txt") %>% as.data.frame()
counts <- fread("/home/yuss/flyG4/result/qians.RNAseq/012.1.counts.txt") %>% as.data.frame()
tpm.calculate = function(exprset,len){
  readperlength = t(do.call(rbind, lapply(1:ncol(exprset), function(i){
    exprset[,i]/len})))
  totalcounts <- colSums(readperlength)
  tpm = t(apply(readperlength, 1, function(x) 10^6 * x/totalcounts)) %>% as.data.frame()
  colnames(tpm) = colnames(exprset)
  row.names(tpm) = row.names(exprset)
  return(tpm)
}
counts$length = gene.length[match(counts$geneid,gene.length$gene.name),2]
counts <- column_to_rownames(counts,var= "geneid")
tpm = tpm.calculate(counts[,-13],counts$length) #ncol() 函数返回矩阵的列数
tpm$DMSO <- rowMeans(tpm[,1:3])
tpm$TMP100 <- rowMeans(tpm[,4:6])
tpm$TMP25 <- rowMeans(tpm[,7:9])
tpm$TMP50 <- rowMeans(tpm[,10:12])
tpm$gene_id <- rownames(tpm)
write.table(tpm,file = paste0("/home/yuss/flyG4/result/qians.RNAseq/",Num,"Tpm.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)

rm(list = ls());gc();rm(list = ls())#清空
Num = "012.1."
#### 3.热图 ####
##相关性热图
library(pheatmap)
tpm <- fread("/home/yuss/flyG4/result/qians.RNAseq/012.1.Tpm.txt") %>% as.data.frame()
rownames(tpm) <- tpm$gene_id
# kc.tpm <- kc.tpm[,c(1:3,7:8)]
tpm <- tpm[,1:12]
# tpm <- as.numeric(tpm) 
table <- cor(tpm, method = "pearson") #cor(x, y, method)用于测量两个向量之间的相关系数值
pheatmap(table, cluster_rows = T, #cluster_rows=T:对行进行集群分析
         cluster_cols = T, #对列进行集群分析
         display_numbers = T,
         number_format = "%.2f") #显示在cell内的数字格式，例如%.2f代表两位小数,默认也是两位小数

data=tpm[which(rowSums(tpm==0)==0),] #每一行都没有元素等于 0 的子集
data <- log2(data+1)
pheatmap(data,  #要绘制热图的矩阵
         color = colorRampPalette(c('blue','white','red'))(100), #热图色块颜色是从蓝到红分为100个等级
         border_color = "black",  #热图中每个色块的边框颜色，NA表示无边框
         scale = "row", #按行进行归一化，"column"表示按列，"none"表示不进行归一化
         cluster_rows = FALSE, #是否对行进行聚类
         cluster_cols = TRUE, #是否对列进行聚类
         legend = TRUE, #是否显示图例
         legend_breaks = c(-1, 0, 1), #设置图例的断点
         legend_labels = c("low","","heigh"), #设置图例断点处的标签
         show_rownames = FALSE, #是否显示行名
         show_colnames = TRUE, #是否显示列名
         fontsize = 8 #字体大小，可以通过fontsize_row、fontsize_col参数分别设置行列名的字体大小
)

