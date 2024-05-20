rm(list = ls());gc()
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
columns(org.Hs.eg.db)  #查看列名
head(keys(org.Hs.eg.db))  #查看行名

#install.packages("readxl")
library("readxl")
genename <- read_excel("/home/yuss/flyG4/30 mass spec proteins.xlsx") 
genename <- genename[,2]
library(stringr)
genename$gene <- str_match(genename$Description, "\\[(.*)\\]")[,2]
genename$gene <- gsub("_HUMAN", "", genename$gene)

##clusterProfiler的 GO 富集
library(clusterProfiler)
#因为表达矩阵是symbol，所以需要转为ENTREZID，才能走clusterProfiler的函数。
gene.df <- bitr(genename$gene, fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = "org.Hs.eg.db")

####GO富集####
ego <- enrichGO(
  gene = gene.df$ENTREZID,
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID",  #指明输入的geneid
  ont = 'ALL',#BP,CC,MF,ALL
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
ego.df <- as.data.frame(ego)
write.table(ego, file = '/home/yuss/flyG4/enrich.go.txt', sep = '\t', row.names = FALSE, quote = FALSE)

barplot(ego)
enrichplot::dotplot(ego)
barplot(ego, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
