#### MLE MSL1####
data <- data.frame(
  Category = c("MLE male", "MLE female", "MSL1 male", "MSL1 female"),
  MergeG4 = c(162,25,114,132),
  shuffle1 = c(149, 7, 137, 70),
  shuffle2 = c(175,4,159,87),
  shuffle3 = c(177,8,177,73),
  shuffle4= c(177,6,150,90),
  shuffle5 =c(182,4,137,70))
data$shuffle1_ratio <- data$MergeG4/data$shuffle1
data$shuffle2_ratio <- data$MergeG4/data$shuffle2
data$shuffle3_ratio <- data$MergeG4/data$shuffle3
data$shuffle4_ratio <- data$MergeG4/data$shuffle4
data$shuffle5_ratio <- data$MergeG4/data$shuffle4

data1 <- data[,c(1,8:12)]
# 使用pivot_longer函数将宽数据转换为长数据，并将第一列作为新的列名
long_data <- pivot_longer(data1, cols = -Category, names_to = "Exam", values_to = "Ratio")
long_data$Category <- factor(long_data$Category,levels=c("MLE female","MLE male","MSL1 female","MSL1 male"))
library(ggplot2)

ggplot(long_data,aes(Category,Ratio))+
  stat_summary(mapping=aes(fill = Category),fun=mean,geom = "bar",fun.args = list(mult=1),width=0.7) + ##柱子
  stat_summary(fun.data=mean_sdl,fun.args = list(mult=1),geom="errorbar",width=0.2)+  # 添加误差线
  labs(x = "",y = "Fold enrichment over random")+
  # stat_compare_means(comparisons = my_comparisons,
  #                  label.y = c(8.5,9),
  #                  aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(0,6)) +
  scale_fill_manual(values = c("#415b9d","#59769e","#7c9d9e","#93b79e")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous(expand = c(0,0))


data <- data.frame(
  Category = c("MLE male", "MLE female", "MSL1 male", "MSL1 female"),
  MergeG4 = c(2620,2757,2680,2663),
  shuffle1 = c(2559,2695,2574,2639),
  shuffle2 = c(2580,2742,2605,2671),
  shuffle3 = c(2472,2635,2485,2577),
  shuffle4= c(2533,2691,2563,2613),
  shuffle5 =c(2494,2661,2534,2602))
data$shuffle1_ratio <- data$MergeG4/data$shuffle1
data$shuffle2_ratio <- data$MergeG4/data$shuffle2
data$shuffle3_ratio <- data$MergeG4/data$shuffle3
data$shuffle4_ratio <- data$MergeG4/data$shuffle4
data$shuffle5_ratio <- data$MergeG4/data$shuffle4
data1 <- data[,c(1,8:12)]
# 使用pivot_longer函数将宽数据转换为长数据，并将第一列作为新的列名
long_data <- pivot_longer(data1, cols = -Category, names_to = "Exam", values_to = "Ratio")
long_data$Category <- factor(long_data$Category,levels=c("MLE female","MLE male","MSL1 female","MSL1 male"))
library(ggplot2)

ggplot(long_data,aes(Category,Ratio))+
  stat_summary(mapping=aes(fill = Category),fun=mean,geom = "bar",fun.args = list(mult=1),width=0.7) + ##柱子
  stat_summary(fun.data=mean_sdl,fun.args = list(mult=1),geom="errorbar",width=0.2)+  # 添加误差线
  labs(x = "",y = "Fold enrichment over random")+
  # stat_compare_means(comparisons = my_comparisons,
  #                  label.y = c(8.5,9),
  #                  aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(0,1.3)) +
  scale_fill_manual(values = c("#415b9d","#59769e","#7c9d9e","#93b79e")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous(expand = c(0,0)) +
  labs(title="Complementary regions of DCC")
