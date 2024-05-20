rm(list = ls());gc();rm(list = ls())#清空
Num="003.5."
#定义输入矩阵
##交集20979
##没有交集43917-20979=22938
##基因组总长度是143726002
##peak长度42473229
##非peak长度是143726002-42473229=101252773
data <- matrix(c(20979, 42473229/143726002*43917, 43917-20979, 101252773/143726002*43917), nrow = 2) 
colnames(data) <- c("peak","not_peak")
rownames(data) <- c("pqs","predict_pqs")
#fisher.test检验
#fisher.test(data) #报错输入值太大
#fisher.test(data, simulate.p.value = TRUE, B = 10000) #用蒙特卡洛（Monte Carlo）算法来计算p值，这种方法通过随机抽样来近似计算 Fisher 检验的 p 值，可以处理更大的输入值，但还是报错
rowSums(data)
#卡方检验
chisq.test(data)

#fisher.test需要数据是整数
##pqs 与 所有peak取交集
data1 <- matrix(c(20979,12978,22938,30939), nrow = 2)
colnames(data1) <- c("Peak regions","Other regions")
rownames(data1) <- c("Observed_PQS","Expected_PQS")
fisher.test(data1)
#数据之前是个matrix，转换数据之前需要转成dataframe
data1_df <- as.data.frame(data1)
#宽数据转为长数据
data1_long <- data1_df %>%
  rownames_to_column(var = "type") %>%
  gather(key = "peak", value = "value", -type)
library(ggplot2)
ggplot(data1_long,aes(x=type,y=value,fill=peak))+
  geom_bar(stat="identity",position="fill",width = 0.6,color='black')+
  theme_classic()+
  xlab("")+
  ylab("")+
  labs(fill="")+
  ggtitle("Fisher's Exact Test , p-value<2.2e-16")+
  coord_flip()+
  scale_x_discrete(labels = c("Observed_PQS" = "Observed\nPQS", "Expected_PQS" = "Expected\nPQS")) +  # 修改 x 坐标轴文本名称
  scale_fill_manual(values = c("#05B9E2","#F27970"))+
  theme(axis.text = element_text (size = 12,color = "black"),
        plot.title = element_text(size = 14),  # 调整标题字体大小
        legend.text = element_text(size = 12))  # 调整图例文本字体大小)

#ggsave("/home/wangdy/BRE_celegans/plot/003_Genage_fisher_test.pdf",width = 4.5,height=3)
ggsave(filename = paste0("/home/yuss/flyG4/result/data.reliability/Picture/",Num,"PeakFisherTest.pdf"),
       device = "pdf",width = 5.5, height = 3) 
##所有peak 与pqs取交集
data1 <- matrix(round(c(60917,42473229/143726002*185156,124239,101252773/143726002*185156),0), nrow = 2)
colnames(data1) <- c("peak","not_peak")
rownames(data1) <- c("pqs","predict_pqs")
fisher.test(data1)



