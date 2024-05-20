# 含有eG4和不含有eG4的基因的表达变化(log2FoldChange)
rm(list = ls());gc();rm(list = ls())
Num = "010.4."
#### kc.Phen ####
kc.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.DE.kc.Phen.txt") %>% as.data.frame()
gene.kc.s2 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.2.gene.kc.s2.txt") %>% as.data.frame()
gene.kc.s2$kc.eG4 <- ifelse(gene.kc.s2$kc==1,"eG4","no eG4")
gene.kc.s2$s2.eG4 <- ifelse(gene.kc.s2$s2==2,"eG4","no eG4")
gene.kc.s2 <- subset(gene.kc.s2,chr=='2L'|chr=='2R'|chr=='3L'|chr=='3R'|chr=='X')
gene.kc.s2$Fc.kcPhen <- kc.Phen[match(gene.kc.s2$id,kc.Phen$geneid),2]
df <- na.omit(gene.kc.s2)

df$chr1 <- ifelse(df$chr=="X","X","A")

library(vioplot)
my_comparisons = list(c("eG4","no eG4"))
ggplot(data = df,aes(x=kc.eG4,y=Fc.kcPhen,fill=kc.eG4)) +
  geom_violin(position = position_dodge(width = 1), scale = 'width') +
  geom_boxplot(position = position_dodge(width = 1), outlier.size = 0.6, width = 0.2, show.legend = FALSE) +
  stat_compare_means(comparisons = my_comparisons,
                   label.y = c(5,5),
                   aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(-3,6)) +
  # scale_fill_manual(values = c("#A4A4A4","#3886B9","#D24E5D")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2](PhenDC3/DMSO))) +
  labs(title="Kc167 cell")

ggplot(data = df,aes(x=kc.eG4,y=Fc.kcPhen,fill=kc.eG4)) +
  geom_violin(position = position_dodge(width = 1), scale = 'width') +
  geom_boxplot(position = position_dodge(width = 1), outlier.size = 0.6, width = 0.2, show.legend = FALSE) +
  facet_grid(~chr1) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(5,5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(-3,6)) +
  # scale_fill_manual(values = c("#A4A4A4","#3886B9","#D24E5D")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2](PhenDC3/DMSO))) +
  labs(title="Kc167 cell")


tapply(df$Fc.kcPhen,df$kc.eG4, mean)
tapply(df$Fc.kcPhen,df$kc.eG4, median)
given_value <- 0
sample_data <- df[df$kc.eG4=="eG4",12]
# 单样本 t 检验(看大于0是否有显著性)
t.test(sample_data, mu = given_value)

aggregate(df$Fc.kcPhen, by = list(df$kc.eG4, df$chr1), FUN = mean)
aggregate(df$Fc.kcPhen, by = list(df$kc.eG4, df$chr1), FUN = median)
# Group.1 Group.2           x
# 1     eG4       A  0.09731426
# 2  no eG4       A -0.05031033
# 3     eG4       X  0.12532028
# 4  no eG4       X -0.04003122


given_value <- 0
sample_data <- df[df$kc.eG4=="eG4"&df$chr1=="A",12]
# 单样本 t 检验
t.test(sample_data, mu = given_value)

given_value <- 0
sample_data <- df[df$kc.eG4=="eG4"&df$chr1=="X",12]
# 单样本 t 检验
t.test(sample_data, mu = given_value)

#### kc.PDS####
kc.PDS <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.DE.kc.PDS.txt") %>% as.data.frame()
gene.kc.s2$Fc.kcPDS <- kc.PDS[match(gene.kc.s2$id,kc.PDS$geneid),2]
df <- na.omit(gene.kc.s2)
df$chr1 <- ifelse(df$chr=="X","X","A")
ggplot(data = df,aes(x=kc.eG4,y=Fc.kcPDS,fill=kc.eG4)) +
  geom_violin(position = position_dodge(width = 1), scale = 'width') +
  geom_boxplot(position = position_dodge(width = 1), outlier.size = 0.6, width = 0.2, show.legend = FALSE) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(2,5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(-3,3)) +
  # scale_fill_manual(values = c("#A4A4A4","#3886B9","#D24E5D")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2](PDS/DMSO))) +
  labs(title="Kc167 cell")
 
ggplot(data = df,aes(x=kc.eG4,y=Fc.kcPDS,fill=kc.eG4)) +
  geom_violin(position = position_dodge(width = 1), scale = 'width') +
  facet_grid(~chr1) +
  geom_boxplot(position = position_dodge(width = 1), outlier.size = 0.6, width = 0.2, show.legend = FALSE) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(2,2),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(-3,3)) +
  # scale_fill_manual(values = c("#A4A4A4","#3886B9","#D24E5D")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2](PDS/DMSO))) +
  labs(title="Kc167 cell")
 
tapply(df$Fc.kcPDS,df$kc.eG4, mean)
aggregate(df$Fc.kcPDS, by = list(df$kc.eG4, df$chr1), FUN = mean)


#### s2.Phen ####
s2.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.DE.s2.Phen.txt") %>% as.data.frame()
gene.kc.s2$Fc.s2Phen <- s2.Phen[match(gene.kc.s2$id,s2.Phen$geneid),2]
df <- na.omit(gene.kc.s2)
df$chr1 <- ifelse(df$chr=="X","X","A")

ggplot(data = df,aes(x=s2.eG4,y=Fc.s2Phen,fill=s2.eG4)) +
  geom_violin(position = position_dodge(width = 1), scale = 'width') +
  geom_boxplot(position = position_dodge(width = 1), outlier.size = 0.6, width = 0.2, show.legend = FALSE) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(9,14),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(-5,10)) +
  # scale_fill_manual(values = c("#A4A4A4","#3886B9","#D24E5D")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2](PhenDC3/DMSO))) +
  labs(title="S2 cell")

ggplot(data = df,aes(x=s2.eG4,y=Fc.s2Phen,fill=s2.eG4)) +
  geom_violin(position = position_dodge(width = 1), scale = 'width') +
  geom_boxplot(position = position_dodge(width = 1), outlier.size = 0.6, width = 0.2, show.legend = FALSE) +
  facet_grid(~chr1) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(9,14),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(-5,10)) +
  # scale_fill_manual(values = c("#A4A4A4","#3886B9","#D24E5D")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2](PhenDC3/DMSO))) +
  labs(title="S2 cell")

tapply(df$Fc.s2Phen,df$s2.eG4, mean)

given_value <- 0
sample_data <- df[df$s2.eG4=="eG4",14]
# 单样本 t 检验(看大于0是否有显著性)
t.test(sample_data, mu = given_value)


aggregate(df$Fc.s2Phen, by = list(df$s2.eG4, df$chr1), FUN = mean)
given_value <- 0
sample_data <- df[df$s2.eG4=="eG4"&df$chr1=="A",14]
# 单样本 t 检验
t_test_result <- t.test(sample_data, mu = given_value)
print(t_test_result)
p-value = 0.008981

sample_data <- df[df$s2.eG4=="eG4"&df$chr1=="X",14]
# 单样本 t 检验
t_test_result <- t.test(sample_data, mu = given_value)
print(t_test_result)
##p-value = 2.775e-11

# 启动子含有eG4和不含有eG4的基因的表达变化(log2FoldChange)
#### kc启动子 ####
rm(list = ls());gc();rm(list = ls())
Num = "010.4."
promoter.kcG4 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter2000.kc.bed") %>% as.data.frame()
kc.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.DE.kc.Phen.txt") %>% as.data.frame()
promoter.kcG4$pro.Fc.kcPhen <- kc.Phen[match(promoter.kcG4$V4,kc.Phen$geneid),2]
promoter.kcG4$kc.eG4 <- ifelse(promoter.kcG4$V6>0,"eG4","no eG4")
promoter.kcG4 <- subset(promoter.kcG4,V1=='2L'|V1=='2R'|V1=='3L'|V1=='3R'|V1=='X')
df <- na.omit(promoter.kcG4)
df$chr1 <- ifelse(df$V1=="X","X","A")

library(vioplot)
my_comparisons = list(c("eG4","no eG4"))
ggplot(data = df,aes(x=kc.eG4,y=pro.Fc.kcPhen,fill=kc.eG4)) +
  geom_violin(position = position_dodge(width = 1), scale = 'width') +
  geom_boxplot(position = position_dodge(width = 1), outlier.size = 0.6, width = 0.2, show.legend = FALSE) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(5,5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(-3,6)) +
  # scale_fill_manual(values = c("#A4A4A4","#3886B9","#D24E5D")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2](PhenDC3/DMSO))) +
  labs(title="Kc167 cell (Promter eG4)")

ggplot(data = df,aes(x=kc.eG4,y=pro.Fc.kcPhen,fill=kc.eG4)) +
  geom_violin(position = position_dodge(width = 1), scale = 'width') +
  geom_boxplot(position = position_dodge(width = 1), outlier.size = 0.6, width = 0.2, show.legend = FALSE) +
  facet_grid(~chr1) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(5,5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(-3,6)) +
  # scale_fill_manual(values = c("#A4A4A4","#3886B9","#D24E5D")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2](PhenDC3/DMSO))) +
  labs(title="Kc167 cell (Promter eG4)")

tapply(df$pro.Fc.kcPhen,df$kc.eG4, mean)
tapply(df$pro.Fc.kcPhen,df$kc.eG4, median)
aggregate(df$pro.Fc.kcPhen, by = list(df$kc.eG4, df$chr1), FUN = mean)
aggregate(df$pro.Fc.kcPhen, by = list(df$kc.eG4, df$chr1), FUN = median)

given_value <- 0
sample_data <- df[df$kc.eG4=="eG4"&df$chr1=="X",7]
# 单样本 t 检验(看大于0是否有显著性)
t.test(sample_data, mu = given_value)
# p-value = 8.003e-09

sample_data <- df[df$kc.eG4=="eG4"&df$chr1=="A",7]
# 单样本 t 检验(看大于0是否有显著性)
t_test_result <- t.test(sample_data, mu = given_value)
print(t_test_result)
# p-value < 2.2e-16
wilcox_test_result <- wilcox.test(sample_data, mu = given_value)
print(wilcox_test_result)

rm(list = ls());gc();rm(list = ls())
Num = "010.4."
#### s2启动子 ####
promoter.s2G4 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter2000.s2.bed") %>% as.data.frame()
s2.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.DE.s2.Phen.txt") %>% as.data.frame()
promoter.s2G4$pro.Fc.s2Phen <- s2.Phen[match(promoter.s2G4$V4,s2.Phen$geneid),2]
promoter.s2G4$s2.eG4 <- ifelse(promoter.s2G4$V6>0,"eG4","no eG4")
promoter.s2G4 <- subset(promoter.s2G4,V1=='2L'|V1=='2R'|V1=='3L'|V1=='3R'|V1=='X')
df <- na.omit(promoter.s2G4)
df$chr1 <- ifelse(df$V1=="X","X","A")
my_comparisons = list(c("eG4","no eG4"))
ggplot(data = df,aes(x=s2.eG4,y=pro.Fc.s2Phen,fill=s2.eG4)) +
  geom_violin(position = position_dodge(width = 1), scale = 'width') +
  geom_boxplot(position = position_dodge(width = 1), outlier.size = 0.6, width = 0.2, show.legend = FALSE) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(8.5,9),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(-5,10)) +
  # scale_fill_manual(values = c("#A4A4A4","#3886B9","#D24E5D")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2](PhenDC3/DMSO))) +
  labs(title="S2 cell (Promter eG4)")

ggplot(data = df,aes(x=s2.eG4,y=pro.Fc.s2Phen,fill=s2.eG4)) +
  geom_violin(position = position_dodge(width = 1), scale = 'width') +
  geom_boxplot(position = position_dodge(width = 1), outlier.size = 0.6, width = 0.2, show.legend = FALSE) +
  facet_grid(~chr1) +
    stat_compare_means(comparisons = my_comparisons,
                     label.y = c(8.5,9),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(-5,10)) +
  # scale_fill_manual(values = c("#A4A4A4","#3886B9","#D24E5D")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2](PhenDC3/DMSO))) +
  labs(title="S2 cell (Promter eG4)")

tapply(df$pro.Fc.s2Phen,df$s2.eG4, mean)
tapply(df$pro.Fc.s2Phen,df$s2.eG4, median)
# eG4    no eG4 
# 0.1172421 0.0476196
aggregate(df$pro.Fc.s2Phen, by = list(df$s2.eG4, df$chr1), FUN = mean)
aggregate(df$pro.Fc.s2Phen, by = list(df$s2.eG4, df$chr1), FUN = median)

given_value <- 0
sample_data <- df[df$s2.eG4=="eG4",7]
# 单样本 t 检验(看大于0是否有显著性)
t.test(sample_data, mu = given_value)
# p-value < 2.2e-16
wilcox.test(sample_data, mu = given_value)

sample_data <- df[df$s2.eG4=="eG4"&df$chr1=="X",7]
# 单样本 t 检验(看大于0是否有显著性)
t.test(sample_data, mu = given_value)
wilcox.test(sample_data, mu = given_value)
# p-value < 2.2e-16

sample_data <- df[df$s2.eG4=="eG4"&df$chr1=="A",7]
# 单样本 t 检验(看大于0是否有显著性)
t.test(sample_data, mu = given_value)
wilcox.test(sample_data, mu = given_value)
# p-value < 2.2e-16
