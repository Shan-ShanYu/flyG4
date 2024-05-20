rm(list = ls());gc();rm(list = ls())#清空
Num = "001.1."
#安装pqsfinder
#BiocManager::install("pqsfinder")
library(pqsfinder)
library(Biostrings)
setwd("/home/yuss/flyG4/script/")
genome <- readDNAStringSet("/home/yuss/flyG4/data/ref/dmel-all-chromosome-r6.19.fasta") 
# chr2L_pqs <- pqsfinder(genome[[1]],overlapping = FALSE, min_score = 50)
# chr2R_pqs <- pqsfinder(genome[[2]],overlapping = FALSE, min_score = 50)
# chr3L_pqs <- pqsfinder(genome[[3]],overlapping = FALSE, min_score = 50)
# chr3R_pqs <- pqsfinder(genome[[4]],overlapping = FALSE, min_score = 50)
# chr4_pqs <- pqsfinder(genome[[5]],overlapping = FALSE, min_score = 50)
# chrX_pqs <- pqsfinder(genome[[6]],overlapping = FALSE, min_score = 50)
# chrY_pqs <- pqsfinder(genome[[7]],overlapping = FALSE, min_score = 50)
##overlapping如果为 true，则将报告所有重叠的 PQS。
##min_score最低PQS分数。默认值 52 显示最佳 Chambers 等人提供的 G4 测序数据的平衡准确性

chr_pqs <- list()
for (i in 1:7) {
  chr_pqs[[i]] <- pqsfinder(genome[[i]],overlapping = FALSE,min_score = 50)
}

a <- c('2L','2R','3L','3R',4,'X','Y')
#提取信息包括染色体、start、end、width、strand、score(结构稳定性)、seq
list <- list()
# list[[1]] <- data.frame(list(chr_pqs[[1]]@ranges,chr_pqs[[1]]@elementMetadata,seq=as.data.frame(chr_pqs[[1]])))
# list[[1]]$chr <- a[1]
for (i in 1:length(chr_pqs)) {
  list[[i]] <- data.frame(list(chr_pqs[[i]]@ranges,chr_pqs[[i]]@elementMetadata))
  list[[i]]$chr <- a[i]
}
df <- do.call('rbind',list)
##合并数据框，rbind根据行进行合并，就是行的叠加，cbind根据列进行合并，即叠加所有列
df <- df[,c(15,1:3,5,4,6:14)]
write.table(df,file = paste0("/home/yuss/flyG4/result/PQS/",Num,"dmel.pqs.txt"), sep = '\t', row.names = F, quote = FALSE)
