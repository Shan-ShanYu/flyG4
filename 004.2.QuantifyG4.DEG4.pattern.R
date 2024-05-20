rm(list = ls());gc();rm(list = ls())
Num = "004.2."

#### merge_G4 ####
#*表达水平-----------------------------------------------------------------------------------
DEG4 <- fread("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/004.1.QuantifyG4.DEG4.txt")
G4bed <- fread("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/004.1.merge.G4wide1K.bed")
G4bed$id <- 1:nrow(G4bed)
DEG4 <- merge(DEG4,G4bed,by.x="id",by.y="id")
DEG4 %<>% dplyr::filter(!V1 %in% c(4,"Y")) ##V1列的元素不在向量 (4, "Y") 中的行将被保留
DEG4$group %<>% gsub("-biased","",.) ##其中的 . 表示当前正在处理的元素
##baseMean是所有样本标准化后的平均表达水平
my_comparisons <- list(c("Female","Male"),c("Male","Unbias"),c("Female","Unbias"))
ggplot(DEG4,aes(group,log10(baseMean),fill=group))+
  geom_boxplot(outlier.colour = "white",notch = TRUE)+
  ylab(bquote(Mean~expression~level~(Log[10]))) + xlab("Sex bias")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12))+
  coord_cartesian(ylim = c(0,5.8))+
  stat_compare_means(comparisons = my_comparisons,method.args = list(alternative="greater"),
                     label.y = c(4.6,5,5.4),tip.length = 0,size = 4)+
  scale_fill_manual(values = c("#f77aae","#82cbff","grey")) +
  guides(fill="none")
ggsave(filename = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/Picture/",Num,"DEG4.ExpressLevel.pdf"),
        device = "pdf",width = 3.5,height = 3)

#*Distribution分布----------------------------------------------------------------------------
df = table(DEG4$group,DEG4$V1) %>% as.matrix() 
df
df <- apply(df, 1, function(x)x/sum(x)*100) %>% as.data.frame()
df$chr = row.names(df)
df %<>% pivot_longer(cols = 1:3) #将df的列1、2、3进行转换，使其从宽格式变为长格式
ggplot(df, aes(chr,value,fill=name))+
  geom_bar(stat="identity",position = "dodge") + #position = "dodge"：这表示将柱状图进行分组显示，使得不同组的柱子并排排列,即绘制分组柱状图
  theme_bw() +
  scale_fill_manual(values = c("#f77aae","#82cbff","grey")) +
  scale_y_continuous(expand=c(0,0),limits = c(0,35))+ #scale_y_continuous() 可以用于修改数据本身的刻度范围，而 coord_cartesian() 则是调整图形的可视化范围，不会改变数据。
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12),
        legend.position = c(0.5, 0.9),
        legend.direction = "horizontal",
        legend.title = element_blank())+
  ylab("Proportion of eG4 (%)") +
  xlab("")
ggsave(filename = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/Picture/",Num,"DEG4.DistributionPlot.pdf"),
        device = "pdf",width = 4,height = 3)    


df = table(DEG4$group,DEG4$V1) %>% as.matrix() 
df
df <- apply(df, 2, function(x)x/sum(x)*100) %>% as.data.frame()
df_long <- df %>%
  rownames_to_column("Gender") %>%
  gather(Chromosome, Percentage, -Gender) %>%
  mutate(Percentage = as.numeric(Percentage))

ggplot(df_long, aes(Chromosome,Percentage,fill=Gender))+
  geom_bar(stat="identity",position = "dodge") + #position = "dodge"：这表示将柱状图进行分组显示，使得不同组的柱子并排排列,即绘制分组柱状图
  theme_bw() +
  scale_fill_manual(values = c("#f77aae","#82cbff","grey")) +
  scale_y_continuous(expand=c(0,0),limits = c(0,100))+ #scale_y_continuous() 可以用于修改数据本身的刻度范围，而 coord_cartesian() 则是调整图形的可视化范围，不会改变数据。
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12),
        legend.position = c(0.5, 0.9),
        legend.direction = "horizontal",
        legend.title = element_blank())+
  ylab("Proportion of eG4 (%)") +
  xlab("")
ggsave(filename = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/Picture/",Num,"DEG4.DistributionPlot.pdf"),
       device = "pdf",width = 4,height = 3)    

DEG4$chr <- ifelse(DEG4$V1=="X","chrX","chrA")
df = table(DEG4$group,DEG4$chr) %>% as.matrix()
df
df <- apply(df, 1, function(x)x/sum(x)*100) %>% as.data.frame()
df$chr <- row.names(df)
df %<>% pivot_longer(cols = 1:3)
df$percent <- paste0(round(df$value,0),"%")
ggplot(df,aes(chr,value,fill=name)) +
  geom_bar(stat = "identity",position = "dodge") +
  theme_bw() +
  scale_fill_manual(values = c("#f77aae","#82cbff","grey")) +
  scale_y_continuous(expand=c(0,0),limits = c(0,110)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12),
        legend.position = c(0.85, 0.7),
        # legend.direction = "horizontal",
        legend.title = element_blank())+
  ylab("Proportion of G4s (%)") +
  xlab("")+
  labs(title="Merge G4s") +
  geom_text(aes(label = percent),position = position_dodge(width = 0.9), vjust=-0.5, size = 3)
ggsave(filename = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/Picture/",Num,"Merge.DEG4.ChA.DistributionPlot.pdf"),
       device = "pdf",width = 4,height = 3)
#*length长度-----------------------------------------------------------------------------------
raw_merge_g4 <- read.table("/home/yuss/flyG4/result/PQS/001.2.merge.bed",header = T) %>% as.data.frame()
DEG4$raw_start <- raw_merge_g4[match(DEG4$V4,raw_merge_g4$id),2]
DEG4$raw_end <- raw_merge_g4[match(DEG4$V4,raw_merge_g4$id),3]
DEG4$raw_length <- DEG4$raw_end-DEG4$raw_start
tapply(abs(DEG4$raw_length), DEG4$group, median)
library(ggpubr)
ggplot(DEG4,aes(group,raw_length,fill=group)) +
  geom_boxplot(outlier.colour = "white",notch = TRUE) +
  scale_fill_manual(values = c("#f77aae","#82cbff","grey")) +
  xlab("Sex bias") +
  ylab("Length (bp)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12)
  )+
  coord_cartesian(ylim = c(10,60)) +
  stat_compare_means(comparisons = my_comparisons,method.args = list(alternative="greater"), #alternative = "greater" 表示比较方法采用的是单侧的替代假设，即检验第一个分组的均值是否大于第二个分组的均值。
                     label.y = c(49,53,57),tip.length = 0,size = 4) +
  guides(fill = "none")
ggsave(filename = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/Picture/",Num,"DEG4.Length.pdf"),
        device = "pdf",width = 3.5,height = 3)

#*structure结构稳定性---------------------------------------------------------------------
tapply(DEG4$V5, DEG4$group, median)
aggregate(V5 ~ group, DEG4, median)
ggplot(DEG4,aes(group,V5,fill=group)) +
  geom_boxplot(outlier.colour = "white",notch = TRUE) +
  scale_fill_manual(values = c("#f77aae","#82cbff","grey")) +
  xlab("Sex bias") +
  ylab("Structural stability") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12)
  )+
  coord_cartesian(ylim = c(50,110)) +
  stat_compare_means(comparisons = my_comparisons,method.args = list(alternative="greater"), #alternative = "greater" 表示比较方法采用的是单侧的替代假设，即检验第一个分组的均值是否大于第二个分组的均值。
                     label.y = c(90,95,100),tip.length = 0,size = 4) +
  guides(fill = "none")
ggsave(filename = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/Picture/",Num,"DEG4.Structure.pdf"),
        device = "pdf",width = 3.5,height = 3)

#*保守性-----------------------------------------------------------------------------------
phastCons <- fread("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/004.2.phastCons.bed") %>% as.data.frame()
phyloP <- fread("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/004.2.phyloP.bed") %>% as.data.frame()
DEG4$phastCons <- phastCons[match(DEG4$V4,phastCons$V1),6]
DEG4$phyloP <- phyloP[match(DEG4$V4,phyloP$V1),6]
tapply(DEG4$phastCons, DEG4$group, median)
tapply(DEG4$phyloP, DEG4$group, median)
ggplot(DEG4,aes(group,phastCons,fill=group)) +
  geom_boxplot(outlier.colour = "white",notch = TRUE) +
  scale_fill_manual(values = c("#f77aae","#82cbff","grey")) +
  xlab("Sex bias") +
  ylab("PhastCons score") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12)
  )+
  coord_cartesian(ylim = c(0,1.25)) +
  stat_compare_means(comparisons = my_comparisons,method.args = list(alternative="greater"), #alternative = "greater" 表示比较方法采用的是单侧的替代假设，即检验第一个分组的均值是否大于第二个分组的均值。
                     label.y = c(1,1.09,1.18),tip.length = 0,size = 4) +
  guides(fill = "none")
ggsave(filename = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/Picture/",Num,"DEG4.phastCons.pdf"),
        device = "pdf",width = 3.5,height = 3)
ggplot(DEG4,aes(group,phyloP,fill=group)) +
  geom_boxplot(outlier.colour = "white",notch = TRUE) +
  scale_fill_manual(values = c("#f77aae","#82cbff","grey")) +
  xlab("Sex bias") +
  ylab("PhyloP score") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12)
  )+
  coord_cartesian(ylim = c(-1,3.5)) +
  stat_compare_means(comparisons = my_comparisons,method.args = list(alternative="greater"), #alternative = "greater" 表示比较方法采用的是单侧的替代假设，即检验第一个分组的均值是否大于第二个分组的均值。
                     label.y = c(2.6,2.88,3.16),tip.length = 0,size = 4) +
  guides(fill = "none")
ggsave(filename = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/Picture/",Num,"DEG4.PhyloP.pdf"),
        device = "pdf",width = 3.5,height = 3)

#*计算TPM-------------------------------------------------------------------------------------
# TPM（每百万转录本数）：这是一种基于归一化的表达量，通常用于比较不同基因在不同样本中的表达水平。
# TPM 考虑了基因长度和测序深度，以百万转录本数为单位表示表达水平。
G4len <- DEG4[,c("V4","raw_length")]
G4len <- column_to_rownames(G4len,var = "V4")
##上面的是原始eG4长度，但是在做转录组时，eG4已经被划分成等长的1000bp
G4len1 <- G4len
G4len1$raw_length <- "1000"
# 将 raw_length 列的值转换为数值类型
G4len1$raw_length <- as.numeric(G4len1$raw_length)
counts <- fread("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/004.1.QuantifyG4.counts.txt") %>% as.data.frame()
counts$id <- 1:nrow(counts)
counts <- counts[DEG4$id,] %>% as.data.frame()
counts$g4id <- DEG4[match(counts$id,DEG4$id),12]
counts <- counts[,c(6,1,2,3,4)]
class(counts$g4id)
##注意：要先把列转化成向量才能作为行名，那么数据框转为向量的中间一步就是转成矩阵
m=as.matrix(counts[,1])
v = as.vector(m)
row.names(counts) <- v
counts <- counts[,-1]
countsTPM <- function(count, efflength){
  RPK <- count/(efflength/1000)
  PMSC_rpk <- colSums(RPK)/1e6
  RPK/PMSC_rpk
}
tpm <- as.data.frame(apply(counts, 2, countsTPM, efflength=G4len1)) ##数据中的每一列应用countsTPM函数，其中efflength为G4len
colnames(tpm) <- colnames(counts)
tpm$g4id <- rownames(tpm)
tpm <- tpm[,c(5,1,2,3,4)]
write.table(tpm,file = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/",Num,"QuantifyG4.tpm.txt"),
                              col.names=T,row.names=F,quote=F)

##TPM表达水平
rm(list = ls())
Num = "004.2."
tpm <- fread("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/004.2.QuantifyG4.tpm.txt") %>% as.data.frame()

tpm1 <- tpm[,-1]
colnames(tpm1) <- c("Female1","Female2","Male1","Male2")
library(pheatmap)
plot <- pheatmap(tpm1,scale = "row",show_rownames=F) ##对行进行标准化
ggsave(plot,filename = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/Picture/",Num,"mergeG4_TPM.HeatPlot.pdf"),
       device = "pdf",width = 4,height = 3)

DEG4 <- fread("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/004.1.QuantifyG4.DEG4.txt") %>% as.data.frame()
G4bed <- fread("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/004.1.merge.G4wide1K.bed") %>% as.data.frame()
G4bed$id <- 1:nrow(G4bed)
DEG4 <- merge(DEG4,G4bed,by.x="id",by.y="id")
DEG4 %<>% dplyr::filter(!V1 %in% c(4,"Y")) ##V1列的元素不在向量 (4, "Y") 中的行将被保留
DEG4$group %<>% gsub("-biased","",.) ##其中的 . 表示当前正在处理的元素
tpm$group <- DEG4[match(tpm$g4id,DEG4$V4),8]
class(tpm$group)
tpm$sum <- rowSums(tpm[,c("female_rep1","female_rep2","male_rep1","male_rep2")])
tpm$mean <- tpm$sum/4
my_comparisons <- list(c("Female","Male"),c("Male","Unbias"),c("Female","Unbias"))
##平均表达水平
ggplot(tpm,aes(group,log10(mean),fill=group)) +
  geom_boxplot(outlier.color = "white",notch = "TRUE") +
  ylab(bquote(Mean~TPM~(Log[10]))) + xlab("Sex bias")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12))+
  coord_cartesian(ylim = c(-0.5,5)) +
  stat_compare_means(comparisons = my_comparisons,method.args = list(alternative="greater"),
                     label.y = c(3.5,4.0,4.5),tip.length = 0,size = 4)+
  scale_fill_manual(values = c("#f77aae","#82cbff","grey")) +
  guides(fill="none")
ggsave(filename = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/Picture/",Num,"DEG4.MeanTPM.pdf"),
        device = "pdf",width = 3.5,height = 3)
 
##最大表达水平
tpm$max <- apply(tpm[,c("female_rep1","female_rep2","male_rep1","male_rep2")],1,max)
ggplot(tpm,aes(group,log10(max),fill=group)) +
  geom_boxplot(outlier.color = "white",notch = "TRUE") +
  ylab(bquote(Max~TPM~(Log[10]))) + xlab("Sex bias")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12))+
  coord_cartesian(ylim = c(-0.5,5)) +
  stat_compare_means(comparisons = my_comparisons,method.args = list(alternative="greater"),
                     label.y = c(3.5,4.0,4.5),tip.length = 0,size = 4)+
  scale_fill_manual(values = c("#f77aae","#82cbff","grey")) +
  guides(fill="none")
ggsave(filename = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/Picture/",Num,"DEG4.MaxTPM.pdf"),
       device = "pdf",width = 3.5,height = 3)

#### kc_all_G4 ####
DEG4 <- fread("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/004.1.Quantify.kcall_G4.DEG4.txt") %>% as.data.frame()
G4bed <- fread("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/004.1.kc_all.G4wide1K.bed") %>% as.data.frame()
G4bed$id <- 1:nrow(G4bed)
DEG4 <- merge(DEG4,G4bed,by.x="id",by.y="id")
DEG4 %<>% dplyr::filter(!V1 %in% c(4,"Y")) ##V1列的元素不在向量 (4, "Y") 中的行将被保留
DEG4$group %<>% gsub("-biased","",.) ##其中的 . 表示当前正在处理的元素
##baseMean是所有样本标准化后的平均表达水平
#*Distribution分布-----------------------------------------------------------------------
df = table(DEG4$group,DEG4$V1) %>% as.matrix() 
df
df <- apply(df, 1, function(x)x/sum(x)*100) %>% as.data.frame()
df$chr = row.names(df)
df %<>% pivot_longer(cols = 1:3) #将df的列1、2、3进行转换，使其从宽格式变为长格式
ggplot(df, aes(chr,value,fill=name))+
  geom_bar(stat="identity",position = "dodge") + #position = "dodge"：这表示将柱状图进行分组显示，使得不同组的柱子并排排列,即绘制分组柱状图
  theme_bw() +
  scale_fill_manual(values = c("#f77aae","#82cbff","grey")) +
  scale_y_continuous(expand=c(0,0),limits = c(0,35))+ #scale_y_continuous() 可以用于修改数据本身的刻度范围，而 coord_cartesian() 则是调整图形的可视化范围，不会改变数据。
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12),
        legend.position = c(0.5, 0.9),
        legend.direction = "horizontal",
        legend.title = element_blank())+
  ylab("Proportion of G4s (%)") +
  xlab("")+
  labs(title="Kc167 G4s") 
ggsave(filename = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/Picture/",Num,"kc_all.DEG4.DistributionPlot.pdf"),
       device = "pdf",width = 4,height = 3)

DEG4$chr <- ifelse(DEG4$V1=="X","chrX","chrA")
df = table(DEG4$group,DEG4$chr) %>% as.matrix()
df
df <- apply(df, 1, function(x)x/sum(x)*100) %>% as.data.frame()
df$chr <- rownames(df)
df %<>% pivot_longer(cols = 1:3)
df$percent <- paste0(round(df$value,0),"%")

ggplot(df,aes(chr,value,fill=name)) +
  geom_bar(stat = "identity",position = "dodge") +
  theme_bw() +
  scale_fill_manual(values = c("#f77aae","#82cbff","grey")) +
  scale_y_continuous(expand=c(0,0),limits = c(0,110)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12),
        legend.position = c(0.85, 0.7),
        # legend.direction = "horizontal",
        legend.title = element_blank())+
  ylab("Proportion of G4s (%)") +
  xlab("")+
  labs(title="Kc167 G4s") +
  geom_text(aes(label = percent), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) #width = 0.9 表示相邻组之间的距离为 0.9 个单元宽度。
  # geom_text(aes(label = percent), vjust = -0.5, size = 5) #vjust参数用于调整标签的垂直位置
ggsave(filename = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/Picture/",Num,"Kc_all.DEG4.ChA.DistributionPlot.pdf"),
       device = "pdf",width = 4,height = 3)

#### s2_all_G4 ####
DEG4 <- fread("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/004.1.Quantify.s2all_G4.DEG4.txt") %>% as.data.frame()
G4bed <- fread("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/004.1.s2_all.G4wide1K.bed") %>% as.data.frame()
G4bed$id <- 1:nrow(G4bed)
DEG4 <- merge(DEG4,G4bed,by.x="id",by.y="id")
DEG4 %<>% dplyr::filter(!V1 %in% c(4,"Y")) ##V1列的元素不在向量 (4, "Y") 中的行将被保留
DEG4$group %<>% gsub("-biased","",.) ##其中的 . 表示当前正在处理的元素
##baseMean是所有样本标准化后的平均表达水平
#*Distribution分布-----------------------------------------------------------------------
df = table(DEG4$group,DEG4$V1) %>% as.matrix() 
df
df <- apply(df, 1, function(x)x/sum(x)*100) %>% as.data.frame()
df$chr = row.names(df)
df %<>% pivot_longer(cols = 1:3) #将df的列1、2、3进行转换，使其从宽格式变为长格式
ggplot(df, aes(chr,value,fill=name))+
  geom_bar(stat="identity",position = "dodge") + #position = "dodge"：这表示将柱状图进行分组显示，使得不同组的柱子并排排列,即绘制分组柱状图
  theme_bw() +
  scale_fill_manual(values = c("#f77aae","#82cbff","grey")) +
  scale_y_continuous(expand=c(0,0),limits = c(0,35))+ #scale_y_continuous() 可以用于修改数据本身的刻度范围，而 coord_cartesian() 则是调整图形的可视化范围，不会改变数据。
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12),
        legend.position = c(0.5, 0.9),
        legend.direction = "horizontal",
        legend.title = element_blank())+
  ylab("Proportion of G4s (%)") +
  xlab("")+
  labs(title="S2 G4s") 
ggsave(filename = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/Picture/",Num,"S2_all.DEG4.DistributionPlot.pdf"),
       device = "pdf",width = 4,height = 3)

DEG4$chr <- ifelse(DEG4$V1=="X","chrX","chrA")
df = table(DEG4$group,DEG4$chr) %>% as.matrix()
df
df <- apply(df, 1, function(x)x/sum(x)*100) %>% as.data.frame()
df$chr <- rownames(df)
df %<>% pivot_longer(cols = 1:3)
df$percent <- paste0(round(df$value,0),"%")
ggplot(df,aes(chr,value,fill=name)) +
  geom_bar(stat = "identity",position = "dodge") +
  theme_bw() +
  scale_fill_manual(values = c("#f77aae","#82cbff","grey")) +
  scale_y_continuous(expand=c(0,0),limits = c(0,110)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12),
        legend.position = c(0.85, 0.7),
        # legend.direction = "horizontal",
        legend.title = element_blank())+
  ylab("Proportion of G4s (%)") +
  xlab("")+
  labs(title="S2 G4s") +
  geom_text(aes(label = percent),position = position_dodge(width = 0.9), vjust=-0.5, size = 3)
ggsave(filename = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/Picture/",Num,"S2_all.DEG4.ChA.DistributionPlot.pdf"),
       device = "pdf",width = 4,height = 3)


rm(list = ls());gc();rm(list = ls())
Num = "004.2."
#### overlap_G4 ####
DEG4 <- fread("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/004.1.Quantify.overlap_G4.DEG4.txt") %>% as.data.frame()
G4bed <- fread("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/004.1.overlap.G4wide1K.bed") %>% as.data.frame()
G4bed$id <- 1:nrow(G4bed)
DEG4 <- merge(DEG4,G4bed,by.x="id",by.y="id")
DEG4 %<>% dplyr::filter(!V1 %in% c(4,"Y")) ##V1列的元素不在向量 (4, "Y") 中的行将被保留
DEG4$group %<>% gsub("-biased","",.) ##其中的 . 表示当前正在处理的元素
##baseMean是所有样本标准化后的平均表达水平
#*Distribution分布-----------------------------------------------------------------------
df = table(DEG4$group,DEG4$V1) %>% as.matrix() 
df
df <- apply(df, 1, function(x)x/sum(x)*100) %>% as.data.frame()
df$chr = row.names(df)
df %<>% pivot_longer(cols = 1:3) #将df的列1、2、3进行转换，使其从宽格式变为长格式
ggplot(df, aes(chr,value,fill=name))+
  geom_bar(stat="identity",position = "dodge") + #position = "dodge"：这表示将柱状图进行分组显示，使得不同组的柱子并排排列,即绘制分组柱状图
  theme_bw() +
  scale_fill_manual(values = c("#f77aae","#82cbff","grey")) +
  scale_y_continuous(expand=c(0,0),limits = c(0,35))+ #scale_y_continuous() 可以用于修改数据本身的刻度范围，而 coord_cartesian() 则是调整图形的可视化范围，不会改变数据。
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12),
        legend.position = c(0.5, 0.9),
        legend.direction = "horizontal",
        legend.title = element_blank())+
  ylab("Proportion of eG4 (%)") +
  xlab("")+
  labs(title="Overlap eG4") 
ggsave(filename = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/Picture/",Num,"overlap.DEG4.DistributionPlot.pdf"),
       device = "pdf",width = 4,height = 3)

df = table(DEG4$group,DEG4$V1) %>% as.matrix() 
df
df <- apply(df, 2, function(x)x/sum(x)*100) %>% as.data.frame()
df_long <- df %>%
  rownames_to_column("Gender") %>%
  gather(Chromosome, Percentage, -Gender) %>%
  mutate(Percentage = as.numeric(Percentage))

ggplot(df_long, aes(Chromosome,Percentage,fill=Gender))+
  geom_bar(stat="identity",position = "dodge") + #position = "dodge"：这表示将柱状图进行分组显示，使得不同组的柱子并排排列,即绘制分组柱状图
  theme_bw() +
  scale_fill_manual(values = c("#f77aae","#82cbff","grey")) +
  scale_y_continuous(expand=c(0,0),limits = c(0,100))+ #scale_y_continuous() 可以用于修改数据本身的刻度范围，而 coord_cartesian() 则是调整图形的可视化范围，不会改变数据。
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12),
        legend.position = c(0.5, 0.9),
        legend.direction = "horizontal",
        legend.title = element_blank())+
  ylab("Proportion of eG4 (%)") +
  xlab("")
# ggsave(filename = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/Picture/",Num,"DEG4.DistributionPlot.pdf"),
       device = "pdf",width = 4,height = 3)    

DEG4$chr <- ifelse(DEG4$V1=="X","chrX","chrA")
df = table(DEG4$group,DEG4$chr) %>% as.matrix()
df
df <- apply(df, 1, function(x)x/sum(x)*100) %>% as.data.frame()
df$chr <- rownames(df)
df %<>% pivot_longer(cols = 1:3)
df$percent <- paste0(round(df$value,0),"%")
ggplot(df,aes(chr,value,fill=name)) +
  geom_bar(stat = "identity",position = "dodge") +
  theme_bw() +
  scale_fill_manual(values = c("#f77aae","#82cbff","grey")) +
  scale_y_continuous(expand=c(0,0),limits = c(0,110)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12),
        legend.position = c(0.85, 0.7),
        # legend.direction = "horizontal",
        legend.title = element_blank())+
  ylab("Proportion of G4s (%)") +
  xlab("")+
  labs(title="Overlap G4s") +
  geom_text(aes(label = percent),position = position_dodge(width = 0.9), vjust=-0.5, size = 3)
ggsave(filename = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/Picture/",Num,"S2_all.DEG4.ChA.DistributionPlot.pdf"),
       device = "pdf",width = 4,height = 3)

#*表达水平-----------------------------------------------------------------------------------
DEG4 <- fread("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/004.1.Quantify.overlap_G4.DEG4.txt")
G4bed <- fread("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/004.1.overlap.G4wide1K.bed")
G4bed$id <- 1:nrow(G4bed)
DEG4 <- merge(DEG4,G4bed,by.x="id",by.y="id")
DEG4 %<>% dplyr::filter(!V1 %in% c(4,"Y")) ##V1列的元素不在向量 (4, "Y") 中的行将被保留
DEG4$group %<>% gsub("-biased","",.) ##其中的 . 表示当前正在处理的元素
##baseMean是所有样本标准化后的平均表达水平
my_comparisons <- list(c("Female","Male"),c("Male","Unbias"),c("Female","Unbias"))
ggplot(DEG4,aes(group,log10(baseMean),fill=group))+
  geom_boxplot(outlier.colour = "white",notch = TRUE)+
  ylab(bquote(Mean~expression~level~(Log[10]))) + xlab("Sex bias")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12))+
  coord_cartesian(ylim = c(0,5.8))+
  stat_compare_means(comparisons = my_comparisons,method.args = list(alternative="greater"),
                     label.y = c(4.6,5,5.4),tip.length = 0,size = 4)+
  scale_fill_manual(values = c("#f77aae","#82cbff","grey")) +
  guides(fill="none")
ggsave(filename = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/Picture/",Num,"overlapDEG4.ExpressLevel.pdf"),
       device = "pdf",width = 3.5,height = 3)

#*length长度-----------------------------------------------------------------------------------
# 感觉没必要计算长度，eG4在匹配表达量时，将eG4前后延伸500bp，再计算之前的eG4长度，觉得不一样
rawg4 <- read.table("/home/yuss/flyG4/result/PQS/001.2.overlap.bed",header = T) %>% as.data.frame()
DEG4$raw_start <- rawg4[match(DEG4$V4,rawg4$id),2]
DEG4$raw_end <- rawg4[match(DEG4$V4,rawg4$id),3]
DEG4$raw_length <- DEG4$raw_end-DEG4$raw_start
tapply(abs(DEG4$raw_length), DEG4$group, median)
library(ggpubr)
ggplot(DEG4,aes(group,raw_length,fill=group)) +
  geom_boxplot(outlier.colour = "white",notch = TRUE) +
  scale_fill_manual(values = c("#f77aae","#82cbff","grey")) +
  xlab("Sex bias") +
  ylab("Length (bp)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12)
  )+
  coord_cartesian(ylim = c(10,60)) +
  stat_compare_means(comparisons = my_comparisons,method.args = list(alternative="greater"), #alternative = "greater" 表示比较方法采用的是单侧的替代假设，即检验第一个分组的均值是否大于第二个分组的均值。
                     label.y = c(49,53,57),tip.length = 0,size = 4) +
  guides(fill = "none")
ggsave(filename = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/Picture/",Num,"overlapDEG4.Length.pdf"),
       device = "pdf",width = 3.5,height = 3)

#*structure结构稳定性---------------------------------------------------------------------
tapply(DEG4$V5, DEG4$group, median)
aggregate(V5 ~ group, DEG4, median)
ggplot(DEG4,aes(group,V5,fill=group)) +
  geom_boxplot(outlier.colour = "white",notch = TRUE) +
  scale_fill_manual(values = c("#f77aae","#82cbff","grey")) +
  xlab("Sex bias") +
  ylab("Structural stability") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12)
  )+
  coord_cartesian(ylim = c(50,110)) +
  stat_compare_means(comparisons = my_comparisons,method.args = list(alternative="greater"), #alternative = "greater" 表示比较方法采用的是单侧的替代假设，即检验第一个分组的均值是否大于第二个分组的均值。
                     label.y = c(90,95,100),tip.length = 0,size = 4) +
  guides(fill = "none")
ggsave(filename = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/Picture/",Num,"overlapDEG4.Structure.pdf"),
       device = "pdf",width = 3.5,height = 3)

#*保守性-----------------------------------------------------------------------------------
phastCons <- fread("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/004.2.overlap.phastCons.bed") %>% as.data.frame()
phyloP <- fread("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/004.2.overlap.phyloP.bed") %>% as.data.frame()
DEG4$phastCons <- phastCons[match(DEG4$V4,phastCons$V1),6]
DEG4$phyloP <- phyloP[match(DEG4$V4,phyloP$V1),6]
tapply(DEG4$phastCons, DEG4$group, median)
tapply(DEG4$phyloP, DEG4$group, median)
ggplot(DEG4,aes(group,phastCons,fill=group)) +
  geom_boxplot(outlier.colour = "white",notch = TRUE) +
  scale_fill_manual(values = c("#f77aae","#82cbff","grey")) +
  xlab("Sex bias") +
  ylab("PhastCons score") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12)
  )+
  coord_cartesian(ylim = c(0,1.3)) +
  stat_compare_means(comparisons = my_comparisons,method.args = list(alternative="greater"), #alternative = "greater" 表示比较方法采用的是单侧的替代假设，即检验第一个分组的均值是否大于第二个分组的均值。
                     label.y = c(1,1.09,1.2),tip.length = 0,size = 4) +
  guides(fill = "none")
ggsave(filename = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/Picture/",Num,"overlapDEG4.phastCons.pdf"),
       device = "pdf",width = 3.5,height = 3)
ggplot(DEG4,aes(group,phyloP,fill=group)) +
  geom_boxplot(outlier.colour = "white",notch = TRUE) +
  scale_fill_manual(values = c("#f77aae","#82cbff","grey")) +
  xlab("Sex bias") +
  ylab("PhyloP score") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10,colour = "black"),axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 10,colour = "black"),axis.title.y = element_text(size = 12)
  )+
  coord_cartesian(ylim = c(-1,4)) +
  stat_compare_means(comparisons = my_comparisons,method.args = list(alternative="greater"), #alternative = "greater" 表示比较方法采用的是单侧的替代假设，即检验第一个分组的均值是否大于第二个分组的均值。
                     label.y = c(2.83,3.2,3.6),tip.length = 0,size = 4) +
  guides(fill = "none")
ggsave(filename = paste0("/home/yuss/flyG4/result/OliverCelniker.Nature.2011.RNAseq/Picture/",Num,"overlapDEG4.PhyloP.pdf"),
       device = "pdf",width = 3.5,height = 3)
