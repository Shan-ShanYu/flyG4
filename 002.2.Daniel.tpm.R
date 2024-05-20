rm(list = ls());gc();rm(list = ls())#清空
Num = "002.2."
#### 非冗余外显子(EXON)长度之和####
#导入数据
exon_length <- fread("/home/yuss/flyG4/data/ref/dmel.nonreExon_length.bed") %>% as.data.frame()
gene.length <- tapply(exon_length$V5,exon_length$V1,sum) %>% as.data.frame()
gene.length1 <- rownames_to_column(gene.length,var = "gene.name")
colnames(gene.length1)[2] <- "length"
write.table(gene.length1,file = paste0("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/",Num,"dmel.genelength.txt"), sep = '\t', col.names = T, row.names = F, quote = FALSE)
counts <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.counts.txt") %>% as.data.frame()
counts <- column_to_rownames(counts,var = "V1")
#tpm
tpm.calculate = function(exprset,len){
  readperlength = t(do.call(rbind, lapply(1:ncol(exprset), function(i){
    exprset[,i]/len})))
  totalcounts <- colSums(readperlength)
  tpm = t(apply(readperlength, 1, function(x) 10^6 * x/totalcounts)) %>% as.data.frame()
  colnames(tpm) = colnames(exprset)
  row.names(tpm) = row.names(exprset)
  return(tpm)
}
counts$length = gene.length[match(row.names(counts),row.names(gene.length)),1]
tpm = tpm.calculate(counts[,-ncol(counts)],counts$length) %>% as.data.frame()#ncol() 函数返回矩阵的列数
colnames(tpm) <- c("kc.tpm","s2.tpm")

#### 分成三类other,eG4,non-eG4 ####
gene.pqs <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.2.gene.pqs.bed") %>% as.data.frame()
colnames(gene.pqs) <- c("chr","start","end","id","strand","gene_symbol","gene.pqs")
gene.pqs$pqs.type <- ifelse(gene.pqs$gene.pqs==0,"other","pqs")
gene.kc.s2 <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.gene.kc.s2.txt") %>% as.data.frame()
gene.pqs$kc <- gene.kc.s2[match(gene.pqs$id,gene.kc.s2$id),7] #匹配kc列
gene.pqs$kc.type <- ifelse(gene.pqs$pqs.type=="pqs",ifelse(gene.pqs$kc=="1","eG4","non-eG4"),"other") #增加kctype列
gene.pqs$s2 <- gene.kc.s2[match(gene.pqs$id,gene.kc.s2$id),8] #匹配s2列
gene.pqs$s2.type <- ifelse(gene.pqs$pqs.type=="pqs",ifelse(gene.pqs$s2=="2","eG4","non-eG4"),"other")
#把tpm匹配到gene.pqs，tpm17494,gene17754,匹配不到的表达量设为0
gene.pqs$kc.tpm <- tpm[match(gene.pqs$id,row.names(tpm)),1]
gene.pqs[is.na(gene.pqs$kc.tpm), 13]=0
gene.pqs$s2.tpm <- tpm[match(gene.pqs$id,row.names(tpm)),2]
gene.pqs[is.na(gene.pqs$s2.tpm), 14]=0

#画图(基因是否含有G4的表达水平)
gene.pqs$kc.type <- factor(gene.pqs$kc.type,levels = c("non-eG4","eG4","other"))
gene.pqs$s2.type <- factor(gene.pqs$s2.type,levels = c("non-eG4","eG4","other"))

my_comparisons = list(c("non-eG4","eG4"),c("eG4","other"))
ggplot(data = gene.pqs,aes(x=kc.type,y=log2(kc.tpm+1),fill=kc.type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(12.5,14),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  coord_cartesian(ylim = c(0,15)) 

ggplot(data = gene.pqs,aes(x=s2.type,y=log2(s2.tpm+1),fill=s2.type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(12.5,14),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#4d4d4d","#D7A49D","#92c5de")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  coord_cartesian(ylim = c(0,15)) 

gene.kc.s2$kc.tpm <- tpm[match(gene.kc.s2$id,row.names(tpm)),1]
gene.kc.s2$s2.tpm <- tpm[match(gene.kc.s2$id,row.names(tpm)),2]
#把tpm匹配到gene.kc.s2，tpm17494,gene17754,匹配不到的表达量设为0
gene.kc.s2[is.na(gene.kc.s2$kc.tpm),10]=0
gene.kc.s2[is.na(gene.kc.s2$s2.tpm),11]=0
gene.kc.s2$kc.type <- ifelse(gene.kc.s2$kc==0,"not contain eG4","eG4")
gene.kc.s2$s2.type <- ifelse(gene.kc.s2$s2==0,"not contain eG4","eG4")
gene.kc.s2$kc.type <- factor(gene.kc.s2$kc.type,levels = c("not contain eG4","eG4"))
gene.kc.s2$s2.type <- factor(gene.kc.s2$s2.type,levels = c("not contain eG4","eG4"))
#画图
my_comparisons = list(c("not contain eG4","eG4"))
ggplot(data = gene.kc.s2,aes(x=kc.type,y=log2(kc.tpm+1),fill=kc.type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = 12.5,
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#4d4d4d","#FDDBC7")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  coord_cartesian(ylim = c(0,15)) 

ggplot(data = gene.kc.s2,aes(x=s2.type,y=log2(s2.tpm+1),fill=s2.type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = 12.5,
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(values = c("#4d4d4d","#D7A49D")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  coord_cartesian(ylim = c(0,15)) 

