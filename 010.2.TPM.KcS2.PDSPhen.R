rm(list = ls());gc();rm(list = ls())#清空
Num = "010.2."
#### 计算TPM ####
#之前计算过基因长度
gene.length <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.2.dmel.genelength.txt") %>% as.data.frame()
counts <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.kc.counts.txt") %>% as.data.frame()
tpm.calculate = function(exprset,len){
  readperlength = t(do.call(rbind, lapply(1:ncol(exprset), function(i){
    exprset[,i]/len})))
  totalcounts <- colSums(readperlength)
  tpm = t(apply(readperlength, 1, function(x) 10^6 * x/totalcounts)) %>% as.data.frame()
  colnames(tpm) = colnames(exprset)
  row.names(tpm) = row.names(exprset)
  return(tpm)
}
counts$length = gene.length[match(counts$geneid,gene.length$gene.name),2]
counts <- column_to_rownames(counts,var= "geneid")
tpm = tpm.calculate(counts[,-9],counts$length) #ncol() 函数返回矩阵的列数
tpm$con <- rowMeans(tpm[,1:3])
tpm$PDS <- rowMeans(tpm[,4:6])
tpm$Phen <- rowMeans(tpm[,7:8])
tpm$gene_id <- rownames(tpm)
write.table(tpm,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"KcTpm.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)

counts <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.s2.counts.txt") %>% as.data.frame()
tpm.calculate = function(exprset,len){
  readperlength = t(do.call(rbind, lapply(1:ncol(exprset), function(i){
    exprset[,i]/len})))
  totalcounts <- colSums(readperlength)
  tpm = t(apply(readperlength, 1, function(x) 10^6 * x/totalcounts)) %>% as.data.frame()
  colnames(tpm) = colnames(exprset)
  row.names(tpm) = row.names(exprset)
  return(tpm)
}
counts$length = gene.length[match(counts$geneid,gene.length$gene.name),2]
counts <- column_to_rownames(counts,var= "geneid")
tpm = tpm.calculate(counts[,-9],counts$length) #ncol() 函数返回矩阵的列数
tpm$con <- rowMeans(tpm[,1:3])
tpm$PDS <- rowMeans(tpm[,4:6])
tpm$Phen <- rowMeans(tpm[,7:8])
tpm$gene_id <- rownames(tpm)
write.table(tpm,file = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/",Num,"S2Tpm.txt"),
            sep = '\t',col.names = T,row.names = F,quote = F)

#### kc细胞含有eG4的基因不同处理下基因表达的比较 ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "010.2."
gene.kc.s2 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.2.gene.kc.s2.txt") %>% as.data.frame()
gene.kc.s2$kc.eG4 <- ifelse(gene.kc.s2$kc==1,"eG4","no eG4")
#把tpm匹配到gene.pqs，tpm17494,gene17754,匹配不到的表达量设为0
tpm <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.2.KcTpm.txt") %>% as.data.frame()
gene.kc.s2$con <- tpm[match(gene.kc.s2$id,tpm$gene_id),9]
gene.kc.s2$PDS <- tpm[match(gene.kc.s2$id,tpm$gene_id),10]
gene.kc.s2$Phen <- tpm[match(gene.kc.s2$id,tpm$gene_id),11]
#剔除
gene.kc.s2 <- subset(gene.kc.s2,kc.eG4=='eG4')

df <- na.omit(gene.kc.s2)
#宽数据转为长数据
long_df <- df %>% gather(key = "treat", value = "kc.tpm", -(1:10))

#画图
long_df$treat <- factor(long_df$treat,levels = c("con","PDS","Phen"))
my_comparisons = list(c("con","PDS"),c("con","Phen"))
ggplot(data = long_df,aes(x=treat,y=log2(kc.tpm+1),fill=treat)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     method="t.test",
                     paired = FALSE,  # 设置为 FALSE，表示独立样本 t 检验
                     label.y = c(12.5,14),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#A4A4A4","#3886B9","#D24E5D")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="eG4 containing genes (Kc167 cell)") +
  coord_cartesian(ylim = c(0,15)) 

long_df$threshold <- ifelse(long_df$kc.tpm>1,"express","no express")
result <- long_df %>%
  group_by(treat) %>%
  count(threshold, name = "count")
wide_df <- spread(result, key = threshold, value = count)
wide_df$sum <- sum(wide_df$express,wide_df$`no express`)
wide_df$ratio <- wide_df$express/wide_df$sum

ggplot(wide_df,aes(x=treat,y=ratio*100,fill=treat))+
  geom_bar(stat = "identity",position = "dodge",width = 0.8)+
  theme_bw()+ 
  ylab("% Expressed genes")+
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=10,colour = "black"),axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,colour = "black"),axis.title.y = element_text(size=12),
        legend.position = "none")+
  scale_y_continuous(expand = c(0,0),limits = c(0,25)) +
  geom_text(aes(label = sprintf("%.0f%%", ratio * 100), vjust = -0.5), color = "black", size = 4)+
  labs(title="eG4 containing genes (Kc167 cell)") +
  scale_fill_manual(values = c("con"="#A4A4A4","PDS"="#3886B9","Phen"="#D24E5D")) +
  scale_x_discrete(labels = c("DMSO","PDS","PhenDC3"))
# ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"Promoter.noeG4.eG4.ExpressRatio_in_KcS2.pdf"),
device = "pdf",width = 4,height = 4.8)

#chrX-------------------------------------------------------
long_df_X <- subset(long_df,chr=='X')
ggplot(data = long_df_X,aes(x=treat,y=log2(kc.tpm+1),fill=treat)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     method="t.test",
                     paired = FALSE,  # 设置为 FALSE，表示独立样本 t 检验
                     label.y = c(12.5,14),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#A4A4A4","#3886B9","#D24E5D")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(#family = "serif", #标题字体
          # face = "bold", #标题加粗
          size = 14),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="eG4 containing genes at X (Kc167)") +
  coord_cartesian(ylim = c(0,15)) +
  scale_x_discrete(labels = c("DMSO","PDS","PhenDC3"))

long_df_X$threshold <- ifelse(long_df_X$kc.tpm>1,"express","no express")
result <- long_df_X %>%
  group_by(treat) %>%
  count(threshold, name = "count")
wide_df <- spread(result, key = threshold, value = count)
wide_df$sum <- sum(wide_df$express,wide_df$`no express`)
wide_df$ratio <- wide_df$express/wide_df$sum

ggplot(wide_df,aes(x=treat,y=ratio*100,fill=treat))+
  geom_bar(stat = "identity",position = "dodge",width = 0.8)+
  theme_bw()+ 
  ylab("% Expressed genes")+
  theme(plot.title = element_text(#family = "serif", #标题字体
                                  # face = "bold", #标题加粗
                                  size = 14),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=10,colour = "black"),axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,colour = "black"),axis.title.y = element_text(size=12),
        legend.position = "none")+
  scale_y_continuous(expand = c(0,0),limits = c(0,25)) +
  geom_text(aes(label = sprintf("%.0f%%", ratio * 100), vjust = -0.5), color = "black", size = 4)+
  labs(title="eG4 containing genes at X (Kc167)") +
  scale_fill_manual(values = c("con"="#A4A4A4","PDS"="#3886B9","Phen"="#D24E5D")) +
  scale_x_discrete(labels = c("DMSO","PDS","PhenDC3"))
# ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"Promoter.noeG4.eG4.ExpressRatio_in_KcS2.pdf"),
device = "pdf",width = 4,height = 4.8)

#### 过滤低表达后 kc细胞含有eG4的基因不同处理下基因表达的比较 ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "010.2."
gene.kc.s2 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.2.gene.kc.s2.txt") %>% as.data.frame()
gene.kc.s2$kc.eG4 <- ifelse(gene.kc.s2$kc==1,"eG4","no eG4")
#把tpm匹配到gene.pqs，tpm17494,gene17754,匹配不到的表达量设为0
tpm <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.2.KcTpm.txt") %>% as.data.frame()
tpm <- tpm[tpm$con>=1&tpm$Phen>=1,]
gene.kc.s2$con <- tpm[match(gene.kc.s2$id,tpm$gene_id),9]
gene.kc.s2$Phen <- tpm[match(gene.kc.s2$id,tpm$gene_id),11]
#剔除
gene.kc.s2 <- subset(gene.kc.s2,kc.eG4=='eG4')
df <- na.omit(gene.kc.s2)
#宽数据转为长数据
long_df <- df %>% gather(key = "treat", value = "kc.tpm", -(1:10))

#画图
long_df$treat <- factor(long_df$treat,levels = c("con","Phen"))
my_comparisons = list(c("con","Phen"))
ggplot(data = long_df,aes(x=treat,y=log2(kc.tpm+1),fill=treat)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(12.5,14),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#A4A4A4","#3886B9","#D24E5D")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="eG4 containing genes (Kc167 cell)") +
  coord_cartesian(ylim = c(0,15)) 

#chrX-------------------------------------------------------
long_df_X <- subset(long_df,chr=='X')
ggplot(data = long_df_X,aes(x=treat,y=log2(kc.tpm+1),fill=treat)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     method="t.test",
                     paired = FALSE,  # 设置为 FALSE，表示独立样本 t 检验
                     label.y = c(12.5,14),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#A4A4A4","#3886B9","#D24E5D")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(#family = "serif", #标题字体
          # face = "bold", #标题加粗
          size = 14),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="eG4 containing genes at X (Kc167)") +
  coord_cartesian(ylim = c(0,15)) +
  scale_x_discrete(labels = c("DMSO","PhenDC3"))

long_df_X$threshold <- ifelse(long_df_X$kc.tpm>1,"express","no express")
result <- long_df_X %>%
  group_by(treat) %>%
  count(threshold, name = "count")
wide_df <- spread(result, key = threshold, value = count)
wide_df$sum <- sum(wide_df$express,wide_df$`no express`)
wide_df$ratio <- wide_df$express/wide_df$sum

ggplot(wide_df,aes(x=treat,y=ratio*100,fill=treat))+
  geom_bar(stat = "identity",position = "dodge",width = 0.8)+
  theme_bw()+ 
  ylab("% Expressed genes")+
  theme(plot.title = element_text(#family = "serif", #标题字体
    # face = "bold", #标题加粗
    size = 14),
    panel.grid.major =element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size=10,colour = "black"),axis.title.x = element_blank(),
    axis.text.y = element_text(size=10,colour = "black"),axis.title.y = element_text(size=12),
    legend.position = "none")+
  scale_y_continuous(expand = c(0,0),limits = c(0,25)) +
  geom_text(aes(label = sprintf("%.0f%%", ratio * 100), vjust = -0.5), color = "black", size = 4)+
  labs(title="eG4 containing genes at X (Kc167)") +
  scale_fill_manual(values = c("con"="#A4A4A4","PDS"="#3886B9","Phen"="#D24E5D")) +
  scale_x_discrete(labels = c("DMSO","PDS","PhenDC3"))
# ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"Promoter.noeG4.eG4.ExpressRatio_in_KcS2.pdf"),
device = "pdf",width = 4,height = 4.8)


#### s2细胞含有eG4的基因不同处理下基因表达的比较 ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "010.2."
gene.kc.s2 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.2.gene.kc.s2.txt") %>% as.data.frame()
gene.kc.s2$s2.eG4 <- ifelse(gene.kc.s2$s2==2,"eG4","no eG4")
#把tpm匹配到gene.pqs，tpm17494,gene17754,匹配不到的表达量设为0
tpm <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.2.S2Tpm.txt") %>% as.data.frame()
gene.kc.s2$con <- tpm[match(gene.kc.s2$id,tpm$gene_id),9]
gene.kc.s2$PDS <- tpm[match(gene.kc.s2$id,tpm$gene_id),10]
gene.kc.s2$Phen <- tpm[match(gene.kc.s2$id,tpm$gene_id),11]
#剔除
gene.kc.s2 <- subset(gene.kc.s2,s2.eG4=='eG4')
# gene.kc.s2 <- subset(gene.kc.s2,chr=='2L'|chr=='2R'|chr=='3L'|chr=='3R'|chr=='X') ##根据eG4剔除后，就只剩下符合的染色体了
df <- na.omit(gene.kc.s2)
#宽数据转为长数据
long_df <- df %>% gather(key = "treat", value = "s2.tpm", -(1:10))

#画图
long_df$treat <- factor(long_df$treat,levels = c("con","PDS","Phen"))
my_comparisons = list(c("con","PDS"),c("con","Phen"))
ggplot(data = long_df,aes(x=treat,y=log2(s2.tpm+1),fill=treat)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     method="t.test",
                     paired = FALSE,  # 设置为 FALSE，表示独立样本 t 检验
                     label.y = c(12.5,14),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#A4A4A4","#3886B9","#D24E5D")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="eG4 containing genes (S2 cell)") +
  coord_cartesian(ylim = c(0,15)) +
  scale_x_discrete(labels = c("DMSO","PDS","PhenDC3"))

long_df$threshold <- ifelse(long_df$s2.tpm>1,"express","no express")
result <- long_df %>%
  group_by(treat) %>%
  count(threshold, name = "count")
wide_df <- spread(result, key = threshold, value = count)
wide_df$sum <- sum(wide_df$express,wide_df$`no express`)
wide_df$ratio <- wide_df$express/wide_df$sum

ggplot(wide_df,aes(x=treat,y=ratio*100,fill=treat))+
  geom_bar(stat = "identity",position = "dodge",width = 0.8)+
  theme_bw()+ 
  ylab("% Expressed genes")+
  theme(plot.title = element_text(#family = "serif", #标题字体
                                  #face = "bold", #标题加粗
                                  size = 14),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=10,colour = "black"),axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,colour = "black"),axis.title.y = element_text(size=12),
        legend.position = "none")+
  scale_y_continuous(expand = c(0,0),limits = c(0,25)) +
  geom_text(aes(label = sprintf("%.0f%%", ratio * 100), vjust = -0.5), color = "black", size = 4)+
  labs(title="eG4 containing genes (S2 cell)") +
  scale_fill_manual(values = c("con"="#A4A4A4","PDS"="#3886B9","Phen"="#D24E5D")) +
  scale_x_discrete(labels = c("DMSO","PDS","PhenDC3"))
# ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"Promoter.noeG4.eG4.ExpressRatio_in_KcS2.pdf"),
       device = "pdf",width = 4,height = 4.8)

#chrX------------------------------------------------------------------------------
long_df_X <- subset(long_df,chr=='X')
ggplot(data = long_df_X,aes(x=treat,y=log2(s2.tpm+1),fill=treat)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     method="t.test",
                     paired = FALSE,  # 设置为 FALSE，表示独立样本 t 检验
                     label.y = c(12.5,14),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#A4A4A4","#3886B9","#D24E5D")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="eG4 containing genes at X (S2)") +
  coord_cartesian(ylim = c(0,15)) +
  scale_x_discrete(labels = c("DMSO","PDS","PhenDC3"))

result <- long_df_X %>%
  group_by(treat) %>%
  count(threshold, name = "count")
wide_df <- spread(result, key = threshold, value = count)
wide_df$sum <- sum(wide_df$express,wide_df$`no express`)
wide_df$ratio <- wide_df$express/wide_df$sum

ggplot(wide_df,aes(x=treat,y=ratio*100,fill=treat))+
  geom_bar(stat = "identity",position = "dodge",width = 0.8)+
  theme_bw()+ 
  ylab("% Expressed genes")+
  theme(plot.title = element_text(#family = "serif", #标题字体
    # face = "bold", #标题加粗
    size = 14),
    panel.grid.major =element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size=10,colour = "black"),axis.title.x = element_blank(),
    axis.text.y = element_text(size=10,colour = "black"),axis.title.y = element_text(size=12),
    legend.position = "none")+
  scale_y_continuous(expand = c(0,0),limits = c(0,25)) +
  geom_text(aes(label = sprintf("%.0f%%", ratio * 100), vjust = -0.5), color = "black", size = 4)+
  labs(title="eG4 containing genes at X (S2 cell)") +
  scale_fill_manual(values = c("con"="#A4A4A4","PDS"="#3886B9","Phen"="#D24E5D")) +
  scale_x_discrete(labels = c("DMSO","PDS","PhenDC3"))
# ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"Promoter.noeG4.eG4.ExpressRatio_in_KcS2.pdf"),
device = "pdf",width = 4,height = 4.8)

rm(list = ls());gc();rm(list = ls())#清空
Num = "010.2."
#### kc细胞含有eG4的基因和不含有eG4的基因的表达水平比较 ####
gene.kc.s2 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.2.gene.kc.s2.txt") %>% as.data.frame()
gene.kc.s2$kc.eG4 <- ifelse(gene.kc.s2$kc==1,"eG4","no eG4")
#把tpm匹配到gene.pqs，tpm17494,gene17754,匹配不到的表达量设为0
tpm <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.2.KcTpm.txt") %>% as.data.frame()
gene.kc.s2$con <- tpm[match(gene.kc.s2$id,tpm$gene_id),9]
gene.kc.s2$PDS <- tpm[match(gene.kc.s2$id,tpm$gene_id),10]
gene.kc.s2$Phen <- tpm[match(gene.kc.s2$id,tpm$gene_id),11]
#剔除
gene.kc.s2 <- subset(gene.kc.s2,chr=='2L'|chr=='2R'|chr=='3L'|chr=='3R'|chr=='X')
df <- na.omit(gene.kc.s2)
df$Phen.con <- log2((df$Phen+1)/(df$con+1))

#画图
# long_df$treat <- factor(long_df$treat,levels = c("con","PDS","Phen"))
my_comparisons = list(c("eG4","no eG4"))
ggplot(data = df,aes(x=kc.eG4,y=df$Phen.con,fill=kc.eG4)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(0.3,14),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(-1,1)) +
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

tapply(df$Phen.con,df$kc.eG4, mean)

df$chr1 <- ifelse(df$chr=="X","X","A")
ggplot(data = df,aes(x=kc.eG4,y=df$Phen.con,fill=kc.eG4)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  facet_grid(~chr1) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(0.3,14),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(-1,1)) +
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


rm(list = ls());gc();rm(list = ls())#清空
Num = "010.2."
#### s2细胞含有eG4的基因和不含有eG4的基因的表达水平比较 ####
gene.kc.s2 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.2.gene.kc.s2.txt") %>% as.data.frame()
gene.kc.s2$s2.eG4 <- ifelse(gene.kc.s2$s2==2,"eG4","no eG4")
#把tpm匹配到gene.pqs，tpm17494,gene17754,匹配不到的表达量设为0
tpm <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.2.S2Tpm.txt") %>% as.data.frame()
gene.kc.s2$con <- tpm[match(gene.kc.s2$id,tpm$gene_id),10]
gene.kc.s2$PDS <- tpm[match(gene.kc.s2$id,tpm$gene_id),11]
gene.kc.s2$Phen <- tpm[match(gene.kc.s2$id,tpm$gene_id),12]
#剔除
gene.kc.s2 <- subset(gene.kc.s2,chr=='2L'|chr=='2R'|chr=='3L'|chr=='3R'|chr=='X')
df <- na.omit(gene.kc.s2)
df$Phen.con <- log2((df$Phen+1)/(df$con+1))

#画图
# long_df$treat <- factor(long_df$treat,levels = c("con","PDS","Phen"))
my_comparisons = list(c("eG4","no eG4"))
ggplot(data = df,aes(x=s2.eG4,y=df$Phen.con,fill=s2.eG4)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(0.8,14),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(-1,2)) +
  # scale_fill_manual(values = c("#A4A4A4","#3886B9","#D24E5D")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2](PDS/DMSO))) +
  labs(title="S2 cell")

tapply(df$Phen.con,df$s2.eG4, mean)

df$chr1 <- ifelse(df$chr=="X","X","A")
ggplot(data = df,aes(x=s2.eG4,y=df$Phen.con,fill=s2.eG4)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  facet_grid(~chr1) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(1,14),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  coord_cartesian(ylim = c(-1,2)) +
  # scale_fill_manual(values = c("#A4A4A4","#3886B9","#D24E5D")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab(bquote(Log[2](PDS/DMSO))) +
  labs(title="S2 cell")

#### 热图 ####
##相关性热图
library(pheatmap)
kc.tpm <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.2.KcTpm.txt") %>% as.data.frame()
s2.tpm <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.2.S2Tpm.txt") %>% as.data.frame()
rownames(kc.tpm) <- kc.tpm$gene_id
kc.tpm <- kc.tpm[,c(1:3,7:8)]
kc.tpm <- as.numeric(kc.tpm)
table <- cor(kc.tpm, method = "pearson") #cor(x, y, method)用于测量两个向量之间的相关系数值
pheatmap(table, cluster_rows = T, #cluster_rows=T:对行进行集群分析
         cluster_cols = T, #对列进行集群分析
         display_numbers = T,
         number_format = "%.2f") #显示在cell内的数字格式，例如%.2f代表两位小数,默认也是两位小数

rownames(s2.tpm) <- s2.tpm$gene_id
s2.tpm <- s2.tpm[,c(1:3,7:8)]
s2.tpm <- as.numeric(s2.tpm)
table <- cor(s2.tpm, method = "pearson") #cor(x, y, method)用于测量两个向量之间的相关系数值
pheatmap(table, cluster_rows = T, #cluster_rows=T:对行进行集群分析
         cluster_cols = T, #对列进行集群分析
         display_numbers = T,
         number_format = "%.2f") #显示在cell内的数字格式，例如%.2f代表两位小数,默认也是两位小数

data=kc.tpm[which(rowSums(kc.tpm==0)==0),] #每一行都没有元素等于 0 的子集
data <- log2(data+1)
pheatmap(data,  #要绘制热图的矩阵
         color = colorRampPalette(c('blue','white','red'))(100), #热图色块颜色是从蓝到红分为100个等级
         border_color = "black",  #热图中每个色块的边框颜色，NA表示无边框
         scale = "row", #按行进行归一化，"column"表示按列，"none"表示不进行归一化
         cluster_rows = FALSE, #是否对行进行聚类
         cluster_cols = TRUE, #是否对列进行聚类
         legend = TRUE, #是否显示图例
         legend_breaks = c(-1, 0, 1), #设置图例的断点
         legend_labels = c("low","","heigh"), #设置图例断点处的标签
         show_rownames = FALSE, #是否显示行名
         show_colnames = TRUE, #是否显示列名
         fontsize = 8 #字体大小，可以通过fontsize_row、fontsize_col参数分别设置行列名的字体大小
)
#### 主成分分析pca ####
kc.tpm <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.2.KcTpm.txt") %>% as.data.frame()
s2.tpm <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.2.S2Tpm.txt") %>% as.data.frame()
rownames(kc.tpm) <- kc.tpm$gene_id
kc.tpm <- kc.tpm[,c(12,1:3,7:8)]
rownames(s2.tpm) <- s2.tpm$gene_id
s2.tpm <- s2.tpm[,c(12,1:3,7:8)]
all <- left_join(kc.tpm,s2.tpm,"gene_id")
rownames(all) <- all$gene_id
all <- all[,-1]
all <- all[apply(all, 1, sum)>0,]
data <- t(all) ##转置
data.pca <- prcomp(data,center = TRUE,scale. = TRUE) ##使用prcomp函数做PCA, 数据要进行标准化，参数都是T
summary(data.pca)

pca1.res <- as.data.frame(data.pca$x)
pca1.res$treat=unlist(lapply(strsplit(rownames(pca1.res),"_"),"[[",2))
pca1.res$cell=unlist(lapply(strsplit(rownames(pca1.res),"_"),"[[",1))
pca1.res
##计算每个主成分对方差的解释度
pca1.var <- data.pca$sdev^2
pca1.var.per <- round(pca1.var/sum(pca1.var)*100, 1)
pca1.res$treat <- factor(pca1.res$treat,
                       levels = c("con", "Phen"),
                       labels = c("Control", "PhenDC3"))
pca1.res$cell <- factor(pca1.res$cell,
                        levels = c("kc", "s2"),
                        labels = c("Kc167", "S2"))
ggplot(pca1.res, aes(x=PC1, y=PC2, color=cell, shape=treat)) +
  geom_point(size=3)+
  theme_bw()+
  xlab(paste("PC1 - ", pca1.var.per[1], "% variance", sep=""))+
  ylab(paste("PC2 - ", pca1.var.per[2], "% variance", sep=""))+
  guides(color=guide_legend(title = "Cell"))+
  guides(shape=guide_legend(title = "Treat"))+
  scale_color_manual(values=c("#FF7F50","#104E8B"))+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

ggsave(filename = paste0("/home/yuss/flyG4/result/KcS2.RNAseq/Picture/",Num,"PCA.pdf"),
       device = "pdf",width = 5.6,height = 3.4)  
