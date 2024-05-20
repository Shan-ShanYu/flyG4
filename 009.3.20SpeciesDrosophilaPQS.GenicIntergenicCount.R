#### PQS数量 ####
rm(list = ls());gc();rm(list = ls())
Num = "009.3."

Count <- read.table("/home/yuss/flyG4/result/20SpeciesDrosophilaPQS/009.3.PQSCountGenicIntergenic.txt.txt",header = F,sep = " ")
setwd("/home/yuss/flyG4/result/20SpeciesDrosophilaPQS")
filename <- list.files(pattern = "\\.bed$")
var_name <- gsub('.pqs.bed', '', filename)
var_name1 <- gsub('009.1.', '', var_name)
colnames(Count) <- c("genic","intergenic")
Count <- Count[-1,]
Count$species <- var_name1
Count_long <- pivot_longer(Count, cols = c("genic","intergenic"), names_to = "type", values_to = "Value")
Count_long$Value <- as.numeric(Count_long$Value)
library(ggplot2)
# 绘制柱状图
ggplot(Count_long, aes(x = species, y = Value/10000, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  ylab(expression("Count (" * x10^4 * ")")) +
  scale_fill_manual(values = c("#7FC97F", "#BDACD3"),
                    breaks = c("genic","intergenic"),
                    labels = c("Genic","Intergenic")) + ##修改图例名称（不是图例标题名称）
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 10)) +
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
    legend.title = element_blank(),
    plot.margin = margin(t = 10,  # 顶部边缘距离
                         r = 5,  # 右边边缘距离
                         b = 5,  # 底部边缘距离
                         l = 5)) # 左边边缘距离
ggsave(filename = paste0("/home/yuss/flyG4/result/20SpeciesDrosophilaPQS/Picture/",Num,"20SpeciesDrosophilaPQS.GenicIntergenicCounts.pdf"),
       device = "pdf",width = 7.5,height = 3.8)

#### PQS密度 ####
##数量除以区域的总长度
setwd("/home/yuss/flyG4/result/20SpeciesDrosophilaPQS/009.3.GCContent")
path <- "/home/yuss/flyG4/result/20SpeciesDrosophilaPQS/009.3.GCContent"
files <- list.files(pattern = "\\.txt$")
filespath <- lapply(files, function(x)paste(path,x,sep = '/'))
data <- list()
data <- lapply(filespath, function(x)fread(x))
b <- gsub(".GC.txt"," ",files) ##gsub("目标字符", "替换字符", 对象)

data2 <- list()
for (i in 1:40) {
  data2[[i]] <- data[[i]][,c(1:3,5,12)]
  data2[[i]]$group <- b[i]
  colnames(data2[[i]]) <- c("chr","start","end","gc","length","group")
}

head(data2[[1]]$gc)
mean(data2[[1]]$gc)

length_list <- list()
for (i in 1:40) {
  length_list[[i]] <- data.frame(group = b[i], length=sum(data2[[i]]$length))
}
df_length <- do.call(rbind,length_list)
df_length$group <- gsub(" ","",df_length$group)

Count_long$type1 <- paste(Count_long$species,Count_long$type,sep=".")
Count_long$length <- df_length[match(Count_long$type1,df_length$group),2]
Count_long$lennormal <- Count_long$Value/Count_long$length

ggplot(Count_long, aes(x = species, y = lennormal*1000, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  # ylab(expression(atop("Number of normalized", paste("PQS ", "(", x10^-4, ")")))) +
  ylab("Density (No./kb)") +
  scale_fill_manual(values = c("#7FC97F", "#BDACD3"),
                    breaks = c("genic","intergenic"),
                    labels = c("Genic","Intergenic")) + ##修改图例名称（不是图例标题名称）
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1)) +
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
    legend.title = element_blank(),
    plot.margin = margin(t = 10,  # 顶部边缘距离
                         r = 5,  # 右边边缘距离
                         b = 5,  # 底部边缘距离
                         l = 5)) # 左边边缘距离
ggsave(filename = paste0("/home/yuss/flyG4/result/20SpeciesDrosophilaPQS/Picture/",Num,"20SpeciesDrosophilaPQS.GenicIntergenic.Density.pdf"),
       device = "pdf",width = 7.5,height = 3.8)


#### GC 含量 ####
df_list <- list()
for (i in 1:40) {
  df_list[[i]] <- data.frame(group = b[i], gc = mean(data2[[i]]$gc))
}
df <- do.call(rbind, df_list)
#as.character(head(df$group)) 这列有空格
df$group <- gsub(" ","",df$group)
df2 <- separate(df,group,into=c("species","type"),sep = "\\.")

ggplot(df2, aes(x = species, y = gc, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  ylab(expression("GC Content")) +
  scale_fill_manual(values = c("#7FC97F", "#BDACD3"),
                    breaks = c("genic","intergenic"),
                    labels = c("Genic","Intergenic")) + ##修改图例名称（不是图例标题名称）
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 0.7)) +
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
    legend.title = element_blank(),
    plot.margin = margin(t = 10,  # 顶部边缘距离
                         r = 5,  # 右边边缘距离
                         b = 5,  # 底部边缘距离
                         l = 5)) # 左边边缘距离
ggsave(filename = paste0("/home/yuss/flyG4/result/20SpeciesDrosophilaPQS/Picture/",Num,"20SpeciesDrosophilaPQS.GenicIntergenic.GCContent.pdf"),
       device = "pdf",width = 7.5,height = 3.8)

#### 标准化后的PQS密度 ####
##数量除以区域的总长度（kb），再除以GC含量
Count_long$gc <- df[match(Count_long$type1,df$group),2]
Count_long$normaldensity <- Count_long$lennormal/Count_long$gc
ggplot(Count_long, aes(x = species, y = normaldensity*1000, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  # ylab(expression(atop("Number of normalized", paste("PQS ", "(", x10^-4, ")")))) +
  ylab("Normalized density (No./kb/GC%)") +
  scale_fill_manual(values = c("#7FC97F", "#BDACD3"),
                    breaks = c("genic","intergenic"),
                    labels = c("Genic","Intergenic")) + ##修改图例名称（不是图例标题名称）
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 2)) +
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
    legend.title = element_blank(),
    plot.margin = margin(t = 12,  # 顶部边缘距离
                         r = 5,  # 右边边缘距离
                         b = 5,  # 底部边缘距离
                         l = 5)) # 左边边缘距离
ggsave(filename = paste0("/home/yuss/flyG4/result/20SpeciesDrosophilaPQS/Picture/",Num,"20SpeciesDrosophilaPQS.GenicIntergenic.NormalizedDensity.pdf"),
       device = "pdf",width = 7.5,height = 3.8)
