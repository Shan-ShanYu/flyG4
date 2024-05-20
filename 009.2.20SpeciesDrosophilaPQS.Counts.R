rm(list = ls());gc();rm(list = ls())
Num = "009.2."
setwd("/home/yuss/flyG4/result/20SpeciesDrosophilaPQS")

# 获取所有BED文件
filename <- list.files(pattern = "\\.bed$")

# 存储文件名称和相应的行数
file_info <- data.frame('Filename' = character(length(filename)), 'Count' = numeric(length(filename)))

# 遍历每个BED文件，获取文件名称和行数
for (i in seq_along(filename)) {
  # 提取文件名称
  var_name <- gsub('.pqs.bed', '', filename[i])
  var_name1 <- gsub('009.1.', '', var_name)
  
  # 读取文件行数
  file_data <- read.table(filename[i], sep = '\t', header = FALSE)
  line_count <- nrow(file_data)
  
  # 存储文件名称和行数
  file_info[i, 'Filename'] <- var_name1
  file_info[i, 'Count'] <- line_count
}

# 获取所有55BED文件
filename <- list.files(pattern = "\\.55bed$")

# 存储文件名称和相应的行数
file55_info <- data.frame('Filename' = character(length(filename)), '55Count' = numeric(length(filename)))

# 遍历每个BED文件，获取文件名称和行数
for (i in seq_along(filename)) {
  # 提取文件名称
  var_name <- gsub('.pqs.55bed', '', filename[i])
  var_name1 <- gsub('009.1.', '', var_name)
  
  # 读取文件行数
  file_data <- read.table(filename[i], sep = '\t', header = FALSE)
  line_count <- nrow(file_data)
  
  # 存储文件名称和行数
  file55_info[i, 'Filename'] <- var_name1
  file55_info[i, '55Count'] <- line_count
}

# 获取所有score为60的BED文件
filename <- list.files(pattern = "\\.60bed$")

# 存储文件名称和相应的行数
file60_info <- data.frame('Filename' = character(length(filename)), '60Count' = numeric(length(filename)))

# 遍历每个BED文件，获取文件名称和行数
for (i in seq_along(filename)) {
  # 提取文件名称
  var_name <- gsub('.pqs.60bed', '', filename[i])
  var_name1 <- gsub('009.1.', '', var_name)
  
  # 读取文件行数
  file_data <- read.table(filename[i], sep = '\t', header = FALSE)
  line_count <- nrow(file_data)
  
  # 存储文件名称和行数
  file60_info[i, 'Filename'] <- var_name1
  file60_info[i, '60Count'] <- line_count
}
file_info$score55 <- file55_info[match(file_info$Filename,file55_info$Filename),3]
file_info$score60 <- file60_info[match(file_info$Filename,file60_info$Filename),3]
colnames(file_info)[2] <- "score50"
file_info_long <- pivot_longer(file_info, cols = starts_with("score"), names_to = "Threshold", values_to = "Value")
file_info_long$type <- gsub("score","",file_info_long$Threshold)

file_info_long$Value10000 <- file_info_long$Value/10000
library(ggplot2)
# 绘制柱状图
ggplot(file_info_long, aes(x = Filename, y = Value10000, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  ylab(expression("Count (" * x10^4 * ")")) +
  scale_fill_manual(values = c("#9DC9E1", "#3182BD", "#08519C"),name = "Cutoff") +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 13)) +
  theme_bw() +
  theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 14),
      axis.text.y = element_text(size = 14, colour = "black"),
      axis.title.y = element_text(size = 16, colour = "black"),
      axis.title.x = element_blank(),
      legend.direction = "horizontal", ##设置图例水平放置
      legend.position = c(0.2,0.88),
      legend.text = element_text(size = 13),  # 调整图例文字大小
      legend.title = element_text(size = 13)  # 调整图例标题大小
  )
ggsave(filename = paste0("/home/yuss/flyG4/result/20SpeciesDrosophilaPQS/Picture/",Num,"20SpeciesDrosophilaPQSCounts.pdf"),
       device = "pdf",width = 7.5,height = 3.8)
