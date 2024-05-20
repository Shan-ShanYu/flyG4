cor_matrix <- read.table("/home/yuss/flyG4/result/data.reliability/003.8.kc.scatterplot.SpearmanCorr.readCounts.tab", header = TRUE, row.names = 1)
# 创建一个与相关系数矩阵相同大小的空数据框，用于存储 p 值
# 创建一个与相关系数矩阵相同大小的空数据框，用于存储 p 值
p_value_matrix <- matrix(0, nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))

# 遍历相关系数矩阵的每一对样本
for (i in 1:nrow(cor_matrix)) {
  for (j in 1:ncol(cor_matrix)) {
    # 跳过对角线的重复计算
    if (i != j) {
      # 提取当前样本对的相关系数
      cor_value <- cor_matrix[i, j]
      
      # 计算 p 值
      p_value <- 2 * (1 - pnorm(abs(cor_value)))
      
      # 将计算得到的 p 值填入 p_value_matrix
      p_value_matrix[i, j] <- p_value
    }
  }
}

# 将 p 值矩阵转换成数据框
p_value_df <- as.data.frame(p_value_matrix, row.names = rownames(cor_matrix), col.names = colnames(cor_matrix))


# 提取相关性矩阵中的相关性值
cor_values <- as.vector(cor_matrix)

# 计算 p 值
p_values <- 2 * (1 - pnorm(abs(cor_values)))

# 将 p 值矩阵转换成与相关性矩阵相同的结构
p_value_matrix <- matrix(p_values, nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))

# 将 p 值矩阵转换成数据框
p_value_df <- as.data.frame(p_value_matrix, row.names = rownames(cor_matrix), col.names = colnames(cor_matrix))

rm(list = ls());gc();rm(list = ls())
Num = "005.1."

kc_top_peak <- fread("/home/yuss/flyG4/result/TopMotif/005.1.KC.top1000.narrowPeak") %>% as.data.frame()
s2_top_peak <- fread("/home/yuss/flyG4/result/TopMotif/005.1.S2.top1000.narrowPeak") %>% as.data.frame()
kc_top_peak <- kc_top_peak[,1:3]
colnames(kc_top_peak) <- c("chr","start","end")
chr.list <- c("2L","2R","3L","3R","4","X","Y")
kc_filter <- filter(kc_top_peak, chr %in% chr.list)
s2_top_peak <- s2_top_peak[,1:3]
colnames(s2_top_peak) <- c("chr","start","end")
s2_filter <- filter(s2_top_peak, chr %in% chr.list)
write.table(kc_filter,file = paste0("/home/yuss/flyG4/result/TopMotif/",Num,"KC.top1000.filter.narrowPeak"),
            sep = '\t',col.names = F,row.names = F,quote = F)
write.table(s2_filter,file = paste0("/home/yuss/flyG4/result/TopMotif/",Num,"S2.top1000.filter.narrowPeak"),
            sep = '\t',col.names = F,row.names = F,quote = F)
