#钱师兄做的TMPyP4 和 我做的两个处理数据，看相同的基因在不同样本中的表达是否相同
rm(list = ls());gc();rm(list = ls())
#### 三个处理的差异基因 ####
readcounts <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.counts.txt") %>% as.data.frame()
qianscounts <- fread("/home/yuss/flyG4/result/qians.RNAseq/012.1.counts.txt") %>% as.data.frame() ##TMP_50_rep2离群样本
qianscounts <- qianscounts[,-11] ##剔除TMP_50_rep2离群样本
merged_df <- merge(readcounts, qianscounts, by = "geneid", all.x = F) #all.x = TRUE表示保留df1中的所有行，并在需要时填充df2中的缺失值
rownames(merged_df) <- merged_df$geneid
merged_df <- merged_df[,-1]
sample_name <- factor(colnames(merged_df))
metadata <- data.frame(sample_name)
metadata$sample_name <- as.character(metadata$sample_name)
metadata$treat <- sapply(strsplit(metadata$sample_name, "_"), function(x) paste(x[1],x[2],sep = "_"))
metadata$treat[17] <- "DMSO"
metadata$treat[18] <- "DMSO"
metadata$treat[19] <- "DMSO"
metadata$treat <- as.factor(metadata$treat) ##构建矩阵前要设置因子水平
##1.构建矩阵
dds <- DESeqDataSetFromMatrix(countData = merged_df,
                              colData = metadata,
                              design = ~ treat)
##2.dds标准化
dds <- DESeq(dds)
##3.获取标准化后的数据
normalized_counts <- counts(dds,normalized=T)
##4.有三个重复取平均
normalized_mean_counts = t(apply(normalized_counts, 1, function(a){tapply(a, metadata$treat, mean)}))

TMP100 <- na.omit(as.data.frame(results(dds, contrast = c("treat", "TMP_100","DMSO"), cooksCutoff = FALSE))) #cooksCutoff = FALSE参数是指离群的基因不会用NA表示，之前没设置参数时，将会把离群值设为NA
TMP100$group <- ifelse(TMP100$pvalue<0.05&abs(TMP100$log2FoldChange)>=0.5,ifelse(TMP100$log2FoldChange>0.5,"Up","Down"),"No-sig")
table(TMP100$group)
# Down No-sig     Up 
# 1085   8436    518
TMP100id <- rownames(TMP100)[TMP100$group != "No-sig"]
TMP25 <- na.omit(as.data.frame(results(dds, contrast = c("treat", "TMP_25","DMSO"), cooksCutoff = FALSE))) #cooksCutoff = FALSE参数是指离群的基因不会用NA表示，之前没设置参数时，将会把离群值设为NA
TMP25$group <- ifelse(TMP25$pvalue<0.05&abs(TMP25$log2FoldChange)>=0.5,ifelse(TMP25$log2FoldChange>0.5,"Up","Down"),"No-sig")
table(TMP25$group)
# Down No-sig     Up 
# 712   8958    599
TMP25id <- rownames(TMP25)[TMP25$group != "No-sig"]
TMP50 <- na.omit(as.data.frame(results(dds, contrast = c("treat", "TMP_50","DMSO"), cooksCutoff = FALSE))) #cooksCutoff = FALSE参数是指离群的基因不会用NA表示，之前没设置参数时，将会把离群值设为NA
TMP50$group <- ifelse(TMP50$pvalue<0.05&abs(TMP50$log2FoldChange)>=0.5,ifelse(TMP50$log2FoldChange>0.5,"Up","Down"),"No-sig")
table(TMP50$group)
# Down No-sig     Up 
# 987   8261    330
TMP50id <- rownames(TMP50)[TMP50$group != "No-sig"]
kcPDS <- na.omit(as.data.frame(results(dds, contrast = c("treat", "kc_PDS","kc_con"), cooksCutoff = FALSE))) #cooksCutoff = FALSE参数是指离群的基因不会用NA表示，之前没设置参数时，将会把离群值设为NA
kcPDS$group <- ifelse(kcPDS$pvalue<0.05&abs(kcPDS$log2FoldChange)>=0.5,ifelse(kcPDS$log2FoldChange>0.5,"Up","Down"),"No-sig")
table(kcPDS$group)
# Down No-sig     Up 
# 249   9543     16 
kcPDSid <- rownames(kcPDS)[kcPDS$group != "No-sig"]
kcPhen <- na.omit(as.data.frame(results(dds, contrast = c("treat", "kc_Phen","kc_con"), cooksCutoff = FALSE))) #cooksCutoff = FALSE参数是指离群的基因不会用NA表示，之前没设置参数时，将会把离群值设为NA
kcPhen$group <- ifelse(kcPhen$pvalue<0.05&abs(kcPhen$log2FoldChange)>=0.5,ifelse(kcPhen$log2FoldChange>0.5,"Up","Down"),"No-sig")
table(kcPhen$group)
# Down No-sig     Up 
# 146   9652    932
kcPhenid <- rownames(kcPhen)[kcPhen$group != "No-sig"]
s2PDS <- na.omit(as.data.frame(results(dds, contrast = c("treat", "s2_PDS","s2_con"), cooksCutoff = FALSE))) #cooksCutoff = FALSE参数是指离群的基因不会用NA表示，之前没设置参数时，将会把离群值设为NA
s2PDS$group <- ifelse(s2PDS$pvalue<0.05&abs(s2PDS$log2FoldChange)>=0.5,ifelse(s2PDS$log2FoldChange>0.5,"Up","Down"),"No-sig")
table(s2PDS$group)
# Down No-sig     Up 
# 102   9188     58 
s2PDSid <- rownames(s2PDS)[s2PDS$group != "No-sig"]
s2Phen <- na.omit(as.data.frame(results(dds, contrast = c("treat", "s2_Phen","s2_con"), cooksCutoff = FALSE))) #cooksCutoff = FALSE参数是指离群的基因不会用NA表示，之前没设置参数时，将会把离群值设为NA
s2Phen$group <- ifelse(s2Phen$pvalue<0.05&abs(s2Phen$log2FoldChange)>=0.5,ifelse(s2Phen$log2FoldChange>0.5,"Up","Down"),"No-sig")
table(s2Phen$group)
# Down No-sig     Up 
# 722   9367   1332 
s2Phenid <- rownames(s2Phen)[s2Phen$group != "No-sig"]
##取差异基因的交集
inter <- Reduce(intersect, list(TMP100id, TMP25id, TMP50id, kcPDSid, kcPhenid, s2PDSid, s2Phenid)) ##无
inter <- Reduce(intersect, list(TMP100id, TMP25id, TMP50id, kcPhenid, s2Phenid))

#### TPM 做样本相关性以及聚类####
##匹配交集基因的TPM
qianstpm <- fread("/home/yuss/flyG4/result/qians.RNAseq/012.1.Tpm.txt") %>% as.data.frame()
filtered_qianstpm <- qianstpm[qianstpm$gene_id %in% inter, c(1:10,12,17)]
kc.tpm <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.2.KcTpm.txt") %>% as.data.frame()
s2.tpm <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.2.S2Tpm.txt") %>% as.data.frame()
filtered_kc.tpm <- kc.tpm[kc.tpm$gene_id %in% inter, c(1:3,7,8,12)]
filtered_s2.tpm <- s2.tpm[s2.tpm$gene_id %in% inter, c(1:3,7,8,12)]
##合并
merged_data <- merge(filtered_qianstpm, filtered_kc.tpm, by = "gene_id")
merged_data <- merge(merged_data, filtered_s2.tpm, by = "gene_id")
##可视化的准备
rownames(merged_data) <- merged_data$gene_id
tpm <- merged_data[,2:22]
# tpm <- as.numeric(tpm) 
table <- cor(tpm, method = "pearson") #cor(x, y, method)用于测量两个向量之间的相关系数值
pheatmap(table, cluster_rows = T, #cluster_rows=T:对行进行集群分析
         cluster_cols = T, #对列进行集群分析
         display_numbers = T,
         number_format = "%.2f") #显示在cell内的数字格式，例如%.2f代表两位小数,默认也是两位小数

data <- log2(tpm+1)
pheatmap(data,  #要绘制热图的矩阵
         color = colorRampPalette(c('blue','white','red'))(100), #热图色块颜色是从蓝到红分为100个等级
         border_color = "black",  #热图中每个色块的边框颜色，NA表示无边框
         scale = "row", #按行进行归一化，"column"表示按列，"none"表示不进行归一化
         cluster_rows = TRUE, #是否对行进行聚类
         cluster_cols = TRUE, #是否对列进行聚类
         legend = TRUE, #是否显示图例
         legend_breaks = c(-1, 0, 1), #设置图例的断点
         legend_labels = c("low","","heigh"), #设置图例断点处的标签
         show_rownames = TRUE, #是否显示行名
         show_colnames = TRUE, #是否显示列名
         fontsize = 8 #字体大小，可以通过fontsize_row、fontsize_col参数分别设置行列名的字体大小
)

