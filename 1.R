#### 基因在各个染色体的基因表达水平 ####
res <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.1.DE.kcvss2.txt") %>% as.data.frame()
gene.kc.s2 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.2.gene.kc.s2.txt") %>% as.data.frame()
tpm <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.Tpm.txt") %>% as.data.frame()
gene.kc.s2$kc.tpm <- tpm[match(gene.kc.s2$id,tpm$gene_id),5]
gene.kc.s2$s2.tpm <- tpm[match(gene.kc.s2$id,tpm$gene_id),6]
gene.kc.s2 <- subset(gene.kc.s2,s2.tpm!='NA')
gene.kc.s2$kc.eG4 <- ifelse(gene.kc.s2$kc==1,"eG4","no eG4")
gene.kc.s2$s2.eG4 <- ifelse(gene.kc.s2$s2==2,"eG4","no eG4")
subset_gene.kc.s2 <- subset(gene.kc.s2, chr %in% c("2L", "2R", "3L","3R","X"))
subset_gene.kc.s2 <- subset_gene.kc.s2 %>%
  mutate(chr_type = ifelse(chr %in% c("2L", "2R", "3L", "3R"), "A", "X"))
my_comparisons = list(c("2L","2R"),c("2R","3L"),c("3L","3R"),c("3R","X"))
my_comparisons <- lapply(my_comparisons,as_label)
ggplot(data = subset_gene.kc.s2,aes(x=chr,y=log2(kc.tpm+1),fill=chr)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = rep(c(13,14),2),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  # scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="All genes (Kc167 cells)") +
  coord_cartesian(ylim = c(0,15)) 
#ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"eG4_TPM_in_Kc.pdf"),
device = "pdf",width = 3.5,height = 3)

ggplot(data = subset_gene.kc.s2,aes(x=chr,y=log2(s2.tpm+1),fill=chr)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = rep(c(13,14),2),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  # scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="All genes (S2 cells)") +
  coord_cartesian(ylim = c(0,15)) 

#### 在各个染色体上含有G4基因的比例 ####
df <- table(subset_gene.kc.s2$kc.eG4, subset_gene.kc.s2$chr) %>% as.data.frame()
df <- spread(df,Var1,Freq)
df$sum <- rowSums(df[,2:3])
df$containing_ratio <- df$eG4/df$sum
kc.s2.containing <- df[,c(1,5)]
df <- table(subset_gene.kc.s2$s2.eG4, subset_gene.kc.s2$chr) %>% as.data.frame()
df <- spread(df,Var1,Freq)
df$sum <- rowSums(df[,2:3])
df$containing_ratio <- df$eG4/df$sum
kc.s2.containing$s2 <- df[,5]
colnames(kc.s2.containing)[2] <- "kc"
kc.s2.containing.long <- pivot_longer(kc.s2.containing, cols = c(kc, s2), names_to = "Measure", values_to = "Value")
kc.s2.containing.long <- kc.s2.containing.long %>% mutate(type = paste(Var2, Measure,sep = "_"))
kc.s2.containing.long$type <- factor(kc.s2.containing.long$type,levels = c("2L_kc","2R_kc","3L_kc","3R_kc","X_kc","2L_s2","2R_s2","3L_s2","3R_s2","X_s2" ))
kc.s2.containing.long$Value <- round(kc.s2.containing.long$Value,2)*100
kc.s2.containing.long$percent <- paste(kc.s2.containing.long$Value, "%",sep = "")
ggplot(kc.s2.containing.long,aes(x=type,y=Value,fill=type))+
  geom_bar(stat = "identity",position = "dodge",width = 0.8)+
  theme_bw()+ 
  ylab("% G4 containing genes")+
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=10,colour = "black"),axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,colour = "black"),axis.title.y = element_text(size=12),
        legend.position = "none")+
  scale_y_continuous(expand = c(0,0),limits = c(0,30)) +
  geom_vline(xintercept = 5.5,color="#595959",size=0.4) +
  geom_text(aes(label=percent,vjust=-0.5),color="black",size=4) + 
  labs(title="All genes (Kc167 and S2 cells)") #+
#颜色  scale_fill_manual(values = c("#8EA1CC","#FC8E62","#8EA1CC","#66C3AA"))# +
#  scale_x_discrete(labels = c("A-Kc","eG4","no eG4","eG4"))
#ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"Promoter.noeG4.eG4.ExpressRatio_in_KcS2.pdf"),
device = "pdf",width = 4,height = 4.8)

ggplot(data = subset_gene.kc.s2,aes(x=chr,y=log2(kc.tpm+1),fill=chr)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = rep(c(13,14),2),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  # scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="All genes (Kc167 cells)") +
  coord_cartesian(ylim = c(0,15)) 
#ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"eG4_TPM_in_Kc.pdf"),
device = "pdf",width = 3.5,height = 3)

#### 基因在常染色体和性染色体的基因表达水平 ####
ggplot(data = subset_gene.kc.s2,aes(x=chr_type,y=log2(kc.tpm+1),fill=chr_type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  # stat_compare_means(comparisons = my_comparisons,
  #                    label.y = rep(c(13,14),2),
  #                    aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  # scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="All genes (Kc167 cells)") +
  coord_cartesian(ylim = c(0,15)) 

ggplot(data = subset_gene.kc.s2,aes(x=chr_type,y=log2(s2.tpm+1),fill=chr_type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  # stat_compare_means(comparisons = my_comparisons,
  #                    label.y = rep(c(13,14),2),
  #                    aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  # scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="All genes (S2 cells)") +
  coord_cartesian(ylim = c(0,15)) 

#### 在常染色体和性染色体上含有G4基因的比例 ####
df <- table(subset_gene.kc.s2$kc.eG4, subset_gene.kc.s2$chr_type) %>% as.data.frame()
df <- spread(df,Var1,Freq)
df$sum <- rowSums(df[,2:3])
df$containing_ratio <- df$eG4/df$sum
kc.s2.containing <- df[,c(1,5)]
df <- table(subset_gene.kc.s2$s2.eG4, subset_gene.kc.s2$chr_type) %>% as.data.frame()
df <- spread(df,Var1,Freq)
df$sum <- rowSums(df[,2:3])
df$containing_ratio <- df$eG4/df$sum
kc.s2.containing$s2 <- df[,5]
colnames(kc.s2.containing)[2] <- "kc"
kc.s2.containing.long <- pivot_longer(kc.s2.containing, cols = c(kc, s2), names_to = "Measure", values_to = "Value")
kc.s2.containing.long <- kc.s2.containing.long %>% mutate(type = paste(Var2, Measure,sep = "_"))
kc.s2.containing.long$type <- factor(kc.s2.containing.long$type,levels = c("A_kc","X_kc","A_s2","X_s2" ))
kc.s2.containing.long$Value <- round(kc.s2.containing.long$Value,2)*100
kc.s2.containing.long$percent <- paste(kc.s2.containing.long$Value, "%",sep = "")
ggplot(kc.s2.containing.long,aes(x=type,y=Value,fill=type))+
  geom_bar(stat = "identity",position = "dodge",width = 0.8)+
  theme_bw()+ 
  ylab("% G4 containing genes")+
  theme(plot.title = element_text(family = "serif", #标题字体
                                  face = "bold", #标题加粗
                                  size = 16),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=10,colour = "black"),axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,colour = "black"),axis.title.y = element_text(size=12),
        legend.position = "none")+
  scale_y_continuous(expand = c(0,0),limits = c(0,30)) +
  geom_vline(xintercept = 5.5,color="#595959",size=0.4) +
  geom_text(aes(label=percent,vjust=-0.5),color="black",size=4) + 
  labs(title="All genes (Kc167 and S2 cells)") #+
#颜色  scale_fill_manual(values = c("#8EA1CC","#FC8E62","#8EA1CC","#66C3AA"))# +
#  scale_x_discrete(labels = c("A-Kc","eG4","no eG4","eG4"))
#ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"Promoter.noeG4.eG4.ExpressRatio_in_KcS2.pdf"),
device = "pdf",width = 4,height = 4.8)

ggplot(data = subset_gene.kc.s2,aes(x=chr,y=log2(kc.tpm+1),fill=chr)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = rep(c(13,14),2),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  # scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="All genes (Kc167 cells)") +
  coord_cartesian(ylim = c(0,15)) 
#ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"eG4_TPM_in_Kc.pdf"),
device = "pdf",width = 3.5,height = 3)

#### 比较kc和s2的含有mergeG4在常染色体的表达，性染色体的表达####
gene.mergeG4 <- subset_gene.kc.s2[subset_gene.kc.s2$sum!=0,]
df1 <- gene.mergeG4[,c(1,10,14)]
df1$type <- paste("kc",gene.mergeG4$chr,sep = "_")
df1$class <- paste("kc",gene.mergeG4$chr_type,sep = "_")
colnames(df1)[2] <- "tpm"
df2 <- gene.mergeG4[,c(1,11,14)]
df2$type <- paste("s2",gene.mergeG4$chr,sep = "_")
df2$class <- paste("s2",gene.mergeG4$chr_type,sep = "_")
colnames(df2)[2] <- "tpm"
merge_df <- rbind(df1,df2)

my_comparisons = list(c("kc_A","s2_A"),c("kc_X","s2_X"))
merge_df$class <- factor(merge_df$class,levels = c("kc_A","s2_A","kc_X","s2_X"))
ggplot(data = merge_df,aes(x=class,y=log2(tpm+1),fill=class)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = rep(c(13,14),2),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  # scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="All genes (Kc167 and S2 cells)") +
  coord_cartesian(ylim = c(0,15)) 
#ggsave(filename = paste0("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/Picture/",Num,"eG4_TPM_in_Kc.pdf"),
device = "pdf",width = 3.5,height = 3)

#### ####
gene.onlykcG4 <- subset_gene.kc.s2[subset_gene.kc.s2$kc!="0"&subset_gene.kc.s2$s2=="0",]
df1 <- gene.onlykcG4[,c(1,10)]
df1$type <- "kc"
colnames(df1)[2] <- "tpm"
df2 <- gene.onlykcG4[,c(1,11)]
df2$type <- "s2"
colnames(df2)[2] <- "tpm"
merge_df <- rbind(df1,df2)

ggplot(data = merge_df,aes(x=type,y=log2(tpm+1),fill=type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = list(c("kc","s2")),
                     label.y = 13,
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  # scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="Only containing Kc167 eG4s genes") + ##kcspecific
  coord_cartesian(ylim = c(0,15)) 

gene.onlys2G4 <- subset_gene.kc.s2[subset_gene.kc.s2$sum=="2",]
df1 <- gene.onlys2G4[,c(1,10)]
df1$type <- "kc"
colnames(df1)[2] <- "tpm"
df2 <- gene.onlys2G4[,c(1,11)]
df2$type <- "s2"
colnames(df2)[2] <- "tpm"
merge_df <- rbind(df1,df2)

ggplot(data = merge_df,aes(x=type,y=log2(tpm+1),fill=type)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  stat_compare_means(comparisons = list(c("kc","s2")),
                     label.y = 13,
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  # scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Expression level (TPM)") +
  labs(title="Only containing S2 eG4s genes") +
  coord_cartesian(ylim = c(0,15)) 
