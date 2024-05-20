##上下调表达的基因上，G4的信号
rm(list = ls());gc();rm(list = ls())#清空
Num = "010.3."

##由于基因的长度会影响G4的数量，基因的长度越长，G4的数量可能会更多，如果要做，应该算密度，G4的数量除以基因总长度
##可以看启动子区域，G4的数量，因为长度规定了一样
##可以把G4为0的剃掉，但是要展示有G4和没有G4的饼图

#### 差异基因分成四类看G4的密度 ####
kc.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.DE.kc.Phen.txt") %>% as.data.frame()
kc.Phen$group1 <- ifelse(kc.Phen$pvalue<0.05&abs(kc.Phen$log2FoldChange)>=0.5,ifelse(kc.Phen$log2FoldChange>0.5,ifelse(kc.Phen$pvalue<0.01&kc.Phen$log2FoldChange>1,"Up strong","Up"),ifelse(kc.Phen$pvalue<0.01&kc.Phen$log2FoldChange<=-1,"Down strong","Down")),"No-sig")
table(kc.Phen$group1)
# Down Down strong      No-sig          Up   Up strong 
# 211          13        9648         511         211 
gene.kc <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.gene.kc.bed") %>% as.data.frame()
kc.Phen$num <- gene.kc[match(kc.Phen$geneid,gene.kc$V4),7]
kc.Phen$chr <- gene.kc[match(kc.Phen$geneid,gene.kc$V4),1]
kc.Phen <- kc.Phen[kc.Phen$chr=="X"|kc.Phen$chr=="2L"|kc.Phen$chr=="2R"|kc.Phen$chr=="3L"|kc.Phen$chr=="3R",]
kc.Phen$chr1 <- ifelse(kc.Phen$chr=="X","X","A")
kc.Phen$start <- gene.kc[match(kc.Phen$geneid,gene.kc$V4),2]
kc.Phen$end <- gene.kc[match(kc.Phen$geneid,gene.kc$V4),3]
kc.Phen$length <- kc.Phen$end-kc.Phen$start
kc.Phen$density <- kc.Phen$num/kc.Phen$length

table(kc.Phen$chr1,kc.Phen$group)
data <- matrix((c(189,31,566,149)),nrow = 2)
fisher.test(data) # 0.03111

kc.Phen$group <- factor(kc.Phen$group,levels = c("No-sig","Down","Up"))
my_comparisons = list(c("No-sig","Down"),c("No-sig","Up"))
ggplot(data = kc.Phen,aes(x=group,y=density*1000,fill=group)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(0.45,0.55),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(name = '',
                     values = c('Up'='#D6604D','No-sig'='grey','Down'='#74ADD1'), #手动设置颜色时调整颜色的因子顺序
  )+
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Number of eG4 (per kb)") +
  coord_cartesian(ylim = c(0,1))

tapply(kc.Phen$density, kc.Phen$group, mean)
# No-sig         Down           Up 
# 7.804497e-05 4.944336e-05 2.039283e-04 
ggplot(data = kc.Phen,aes(x=group,y=density*1000,fill=group)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  facet_grid(~chr1) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(0.45,0.55),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(name = '',
                    values = c('Up'='#D6604D','No-sig'='grey','Down'='#74ADD1'), #手动设置颜色时调整颜色的因子顺序
  )+
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Number of eG4 (per kb)") +
  coord_cartesian(ylim = c(0,1))

kc.Phen$group1 <- factor(kc.Phen$group1,levels = c("No-sig","Down","Down strong","Up","Up strong"))
my_comparisons = list(c("No-sig","Down"),c("No-sig","Down strong"),c("No-sig","Up"),c("No-sig","Up strong"))

ggplot(data = kc.Phen,aes(x=group1,y=density*1000,fill=group1)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(0.4,0.45,0.5,0.55),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(name = '',
                    values = c('Up'='#D6604D','No-sig'='grey','Down'='#74ADD1','Up strong'='#ae3b28','Down strong'='#3e8abb'), #手动设置颜色时调整颜色的因子顺序
  )+
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Number of eG4 (per kb)") +
  coord_cartesian(ylim = c(0,1)) 

ggplot(data = kc.Phen,aes(x=group1,y=density*1000,fill=group1)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  facet_grid(~chr1) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(0.4,0.45,0.5,0.55),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(name = '',
                    values = c('Up'='#D6604D','No-sig'='grey','Down'='#74ADD1','Up strong'='#ae3b28','Down strong'='#3e8abb'), #手动设置颜色时调整颜色的因子顺序
  )+
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Number of eG4 (per kb)") +
  coord_cartesian(ylim = c(0,1)) 
table(kc.Phen$chr1, kc.Phen$group1)

kc.Phen$numtype <- ifelse(kc.Phen$num > 0,"have","no")
a <- as.data.frame(table(kc.Phen$group, kc.Phen$numtype))
wide_a <- spread(a, Var2, Freq)
wide_a$sum <- wide_a$have + wide_a$no
wide_a$have.ratio <- wide_a$have/wide_a$sum
wide_a$no.ratio <- wide_a$no/wide_a$sum

ggplot(wide_a, aes(x=Var1,y=have.ratio*100,fill=Var1)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(0.7),width = 0.7,color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值
  coord_cartesian(ylim = c(0,60)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_npg()+ ##更改颜色
  # geom_text(aes(label=have.ratio*100,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题(左下边框，没有网格线)
  ylab("(%) Gene with eG4") + ##通过bquote函数给图标签添加上下标
  geom_signif(y_position=26.5, xmin=1, xmax=2,
              annotation=c("8.8e-04"),tip_length=0)+
  geom_signif(y_position=49, xmin=2, xmax=3,
              annotation=c("2.2e-16"),tip_length=0)+
  geom_signif(y_position=56, xmin=1, xmax=3,
              annotation=c("2.2e-16"),tip_length=0)+
  theme(plot.title = element_text(size = 14),
        axis.title.y = element_text(size = 16), ##y坐标轴标题字体大小
        axis.title.x = element_blank(), ##删除x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.title = element_blank()) +
  labs(title="Kc167 cell")

data <- matrix((c(2242,7312,31,189)),nrow = 2)
fisher.test(data) # p-value = 0.0008801

data <- matrix((c(2242,7312,329,386)),nrow = 2)
fisher.test(data) # p-value < 2.2e-16

data <- matrix((c(31,189,329,386)),nrow = 2)
fisher.test(data) # p-value < 2.2e-16

kc.Phen$numtype <- ifelse(kc.Phen$num > 0,"have","no")
a <- as.data.frame(table(kc.Phen$group1, kc.Phen$numtype))
wide_a <- spread(a, Var2, Freq)
wide_a$sum <- wide_a$have + wide_a$no
wide_a$have.ratio <- wide_a$have/wide_a$sum
wide_a$no.ratio <- wide_a$no/wide_a$sum

ggplot(wide_a, aes(x=Var1,y=have.ratio*100,fill=Var1)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(0.7),width = 0.7,color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值
  coord_cartesian(ylim = c(0,60)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_npg()+ ##更改颜色
  # geom_text(aes(label=have.ratio*100,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题(左下边框，没有网格线)
  ylab("(%) Gene with eG4") + ##通过bquote函数给图标签添加上下标
  geom_signif(y_position=50, xmin=3, xmax=4,
              annotation=c("***"),tip_length=0)+
  geom_signif(y_position=55, xmin=3, xmax=5,
              annotation=c("****"),tip_length=0)+
  theme(plot.title = element_text(size = 14),
        axis.title.y = element_text(size = 16), ##y坐标轴标题字体大小
        axis.title.x = element_blank(), ##删除x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.title = element_blank()) +
  labs(title="Kc167 cell")
##group1
data <- matrix((c(71, 79,314,698)),nrow = 2)
fisher.test(data) # p-value = 0.0001256

data <- matrix((c(71, 79,338,791)),nrow = 2)
fisher.test(data) # p-value = 3.566e-05


#chrX,chA
a <- as.data.frame(table(kc.Phen$group,kc.Phen$chr1, kc.Phen$numtype))
wide_a <- spread(a, Var3, Freq)
wide_a$sum <- wide_a$have + wide_a$no
wide_a$have.ratio <- wide_a$have/wide_a$sum
wide_a$no.ratio <- wide_a$no/wide_a$sum

ggplot(wide_a, aes(x=Var1,y=have.ratio*100,fill=Var2)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(0.7),width = 0.7,color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值
  coord_cartesian(ylim = c(0,60)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_npg()+ ##更改颜色
  # geom_text(aes(label=have.ratio*100,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题(左下边框，没有网格线)
  ylab("(%) Gene with eG4") + ##通过bquote函数给图标签添加上下标
  geom_signif(y_position=57.5, xmin=1.2, xmax=3.2,
              annotation=c("4.5e-10"),tip_length=0)+ #p-value = 0.4004
  geom_signif(y_position=28, xmin=1.2, xmax=2.2,
              annotation=c("0.2992"),tip_length=0)+ #p-value = 0.4004
  geom_signif(y_position=52.5, xmin=2.2, xmax=3.2,
              annotation=c("3.1e-04"),tip_length=0)+ #p-value = 0.4004
  theme(plot.title = element_text(size = 14),
        axis.title.y = element_text(size = 16), ##y坐标轴标题字体大小
        axis.title.x = element_blank(), ##删除x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.title = element_blank()) +
  labs(title="Kc167 cell")

data <- matrix((c(1822,6104,420,1208)),nrow = 2)
fisher.test(data) # 0.01602
data <- matrix((c(26,163,5,26)),nrow = 2)
fisher.test(data) # 0.7804
data <- matrix((c(253,313,76,73)),nrow = 2)
fisher.test(data) # 0.1959

data <- matrix((c(420,1208,76,73)),nrow = 2)
fisher.test(data) # 4.469e-10
data <- matrix((c(420,1208,5,26)),nrow = 2)
fisher.test(data) # 0.2992

data <- matrix((c(76,73,5,26)),nrow = 2)
fisher.test(data) # 0.0003125

##group1  x染色体上只有一个，这样划分不是很准确了
a <- as.data.frame(table(kc.Phen$group1,kc.Phen$chr1, kc.Phen$numtype))
wide_a <- spread(a, Var3, Freq)
wide_a$sum <- wide_a$have + wide_a$no
wide_a$have.ratio <- wide_a$have/wide_a$sum
wide_a$no.ratio <- wide_a$no/wide_a$sum


#### 差异基因分成四类看启动子G4的数量 ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "010.3."
promoter.kcG4 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter2000.kc.bed") %>% as.data.frame()
kc.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.DE.kc.Phen.txt") %>% as.data.frame()
kc.Phen$group1 <- ifelse(kc.Phen$pvalue<0.05&abs(kc.Phen$log2FoldChange)>=0.5,ifelse(kc.Phen$log2FoldChange>0.5,ifelse(kc.Phen$pvalue<0.01&kc.Phen$log2FoldChange>1,"Up strong","Up"),ifelse(kc.Phen$pvalue<0.01&kc.Phen$log2FoldChange<=-1,"Down strong","Down")),"No-sig")
table(kc.Phen$group1)
kc.Phen$promoter.num <- promoter.kcG4[match(kc.Phen$geneid,promoter.kcG4$V4),6]
kc.Phen$chr <- promoter.kcG4[match(kc.Phen$geneid,promoter.kcG4$V4),1]
kc.Phen <- kc.Phen[kc.Phen$chr=="X"|kc.Phen$chr=="2L"|kc.Phen$chr=="2R"|kc.Phen$chr=="3L"|kc.Phen$chr=="3R",]
kc.Phen$chr1 <- ifelse(kc.Phen$chr=="X","X","A")
kc.Phen$group <- factor(kc.Phen$group,levels = c("No-sig","Down","Up"))
my_comparisons = list(c("No-sig","Down"),c("No-sig","Up"))
ggplot(data = kc.Phen,aes(x=group,y=promoter.num,fill=group)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(3,3.5,4,4.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  # scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Number of Promoter eG4") +
  coord_cartesian(ylim = c(0,5)) +
  scale_y_continuous(expand = c(0,0)) 

tapply(kc.Phen$promoter.num, kc.Phen$group, mean)

ggplot(data = kc.Phen,aes(x=group,y=promoter.num,fill=group)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  facet_grid(~chr1)+
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(3,3.5,4,4.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  # scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Number of Promoter eG4") +
  coord_cartesian(ylim = c(0,5)) +
  scale_y_continuous(expand = c(0,0)) 

kc.Phen$group1 <- factor(kc.Phen$group1,levels = c("No-sig","Down","Down strong","Up","Up strong"))
my_comparisons = list(c("No-sig","Down"),c("No-sig","Down strong"),c("No-sig","Up"),c("No-sig","Up strong"))
ggplot(data = kc.Phen,aes(x=group1,y=promoter.num,fill=group1)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(3,3.5,4,4.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  # scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  cowplot::theme_half_open()+
  theme(axis.text.x = element_text(size = 14,angle = 90, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Number of eG4") +
  coord_cartesian(ylim = c(0,5))# +
  # scale_y_continuous(expand = c(0,0)) 

ggplot(data = kc.Phen,aes(x=group1,y=promoter.num,fill=group1)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  facet_grid(~chr1)+
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(3,3.5,4,4.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  # scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  cowplot::theme_half_open()+
  theme(axis.text.x = element_text(size = 14,angle = 90, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Number of eG4") +
  coord_cartesian(ylim = c(0,5))# +
scale_y_continuous(expand = c(0,0)) 

##Down 缺失G4
aggregate(kc.Phen$promoter.num, by = list(kc.Phen$group1, kc.Phen$chr1), FUN = mean)
##把0过滤
kc.Phen.filter <- kc.Phen[kc.Phen$promoter.num != 0, ]
kc.Phen.filter$group1 <- factor(kc.Phen.filter$group1,levels = c("No-sig","Down","Down strong","Up","Up strong"))
my_comparisons = list(c("No-sig","Down"),c("No-sig","Down strong"),c("No-sig","Up"),c("No-sig","Up strong"))
ggplot(data = kc.Phen.filter,aes(x=group1,y=promoter.num,fill=group1)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(3,3.5,4,4.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(name = '',
                    values = c('Up'='#D6604D','No-sig'='grey','Down'='#74ADD1','Up strong'='#ae3b28','Down strong'='#3e8abb') )+#手动设置颜色时调整颜色的因子顺序
  cowplot::theme_half_open()+
  theme(axis.text.x = element_text(size = 14,angle = 90, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Number of eG4 in Promoter") +
  coord_cartesian(ylim = c(0,5)) +
  scale_y_continuous(expand = c(0,0)) 

tapply(kc.Phen.filter$promoter.num, kc.Phen.filter$group1, mean)
# No-sig        Down Down strong          Up   Up strong 
# 1.306122    1.105263    2.000000    1.533708    1.594595 

ggplot(data = kc.Phen.filter,aes(x=group1,y=promoter.num,fill=group1)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  facet_grid(~chr1) +
  # stat_compare_means(comparisons = my_comparisons,
  #                    label.y = c(3,3.5,4,4.5),
  #                    aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(name = '',
                    values = c('Up'='#D6604D','No-sig'='grey','Down'='#74ADD1','Up strong'='#ae3b28','Down strong'='#3e8abb') )+#手动设置颜色时调整颜色的因子顺序
  cowplot::theme_half_open()+
  theme(axis.text.x = element_text(size = 14,angle = 90, hjust = 0.5),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Number of eG4 in Promoter") +
  coord_cartesian(ylim = c(0,5)) +
  scale_y_continuous(expand = c(0,0)) 

tapply(kc.Phen.filter$promoter.num, kc.Phen.filter$group1, mean)

kc.Phen.filter$group <- factor(kc.Phen.filter$group,levels = c("No-sig","Down","Up"))
my_comparisons = list(c("No-sig","Down"),c("No-sig","Up"))
ggplot(data = kc.Phen.filter,aes(x=group,y=promoter.num,fill=group)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(3,3.5,4,4.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(name = '',
                    values = c('Up'='#D6604D','No-sig'='grey','Down'='#74ADD1'), #手动设置颜色时调整颜色的因子顺序
  )+
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Number of eG4 in Promoter") +
  coord_cartesian(ylim = c(0,5)) +
  scale_y_continuous(expand = c(0,0))

tapply(kc.Phen.filter$promoter.num, kc.Phen.filter$group, mean)

ggplot(data = kc.Phen.filter,aes(x=group,y=promoter.num,fill=group)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  facet_grid(~chr1) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(3,3.5,4,4.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(name = '',
                    values = c('Up'='#D6604D','No-sig'='grey','Down'='#74ADD1'), #手动设置颜色时调整颜色的因子顺序
  )+
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Number of eG4 in Promoter") +
  coord_cartesian(ylim = c(0,5)) +
  scale_y_continuous(expand = c(0,0))

kc.Phen$promter.numtype <- ifelse(kc.Phen$promoter.num > 0,"have","no")
a <- as.data.frame(table(kc.Phen$group, kc.Phen$promter.numtype))
wide_a <- spread(a, Var2, Freq)
wide_a$sum <- wide_a$have + wide_a$no
wide_a$have.ratio <- wide_a$have/wide_a$sum
wide_a$no.ratio <- wide_a$no/wide_a$sum

ggplot(wide_a, aes(x=Var1,y=have.ratio*100,fill=Var1)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(0.7),width = 0.7,color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值
  coord_cartesian(ylim = c(0,40)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_npg()+ ##更改颜色
  # geom_text(aes(label=have.ratio*100,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题(左下边框，没有网格线)
  ylab("(%) Gene with Promoter eG4") + ##通过bquote函数给图标签添加上下标
  geom_signif(y_position=26, xmin=1, xmax=2,
              annotation=c("0.03068"),tip_length=0)+
  geom_signif(y_position=36, xmin=2, xmax=3,
              annotation=c("4.9e-07"),tip_length=0)+
  geom_signif(y_position=38, xmin=1, xmax=3,
              annotation=c("1.5e-10"),tip_length=0)+
  theme(plot.title = element_text(size = 14),
        axis.title.y = element_text(size = 16), ##y坐标轴标题字体大小
        axis.title.x = element_blank(), ##删除x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.title = element_blank()) +
  labs(title="Kc167 cell")

data <- matrix((c(2303,7251,39,181)),nrow = 2)
fisher.test(data) # 0.03068
data <- matrix((c(2303,7251,252,463)),nrow = 2)
fisher.test(data) # 1.531e-10
data <- matrix((c(252,463,39,181)),nrow = 2)
fisher.test(data) # 4.934e-07

#chrX,chA
a <- as.data.frame(table(kc.Phen$group,kc.Phen$chr1, kc.Phen$promter.numtype))
wide_a <- spread(a, Var3, Freq)
wide_a$sum <- wide_a$have + wide_a$no
wide_a$have.ratio <- wide_a$have/wide_a$sum
wide_a$no.ratio <- wide_a$no/wide_a$sum

ggplot(wide_a, aes(x=Var1,y=have.ratio*100,fill=Var2)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(0.7),width = 0.7,color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值
  coord_cartesian(ylim = c(0,40)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_npg()+ ##更改颜色
  # geom_text(aes(label=have.ratio*100,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题(左下边框，没有网格线)
  ylab("(%) Gene with Promoter eG4") + ##通过bquote函数给图标签添加上下标
  # geom_signif(y_position=28, xmin=2.2, xmax=3.2,
  #             annotation=c("****"),tip_length=0)+ #p-value = 0.4004
  theme(plot.title = element_text(size = 14),
        axis.title.y = element_text(size = 16), ##y坐标轴标题字体大小
        axis.title.x = element_blank(), ##删除x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.title = element_blank()) +
  labs(title="Kc167 cell")

data <- matrix((c(416,1212,9,22)),nrow = 2)
fisher.test(data) # p-value = 0.6788
data <- matrix((c(416,1212,43,106)),nrow = 2)
fisher.test(data) # p-value = 0.3797
data <- matrix((c(43,106,9,22)),nrow = 2)
fisher.test(data) # p-value = 1


rm(list = ls());gc();rm(list = ls())#清空
Num = "010.3."
#### 差异基因分成四类看G4的密度 ####
s2.phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.DE.s2.Phen.txt") %>% as.data.frame()
s2.phen$group1 <- ifelse(s2.phen$pvalue<0.05&abs(s2.phen$log2FoldChange)>=0.5,ifelse(s2.phen$log2FoldChange>0.5,ifelse(s2.phen$pvalue<0.01&s2.phen$log2FoldChange>1,"Up strong","Up"),ifelse(s2.phen$pvalue<0.01&s2.phen$log2FoldChange<=-1,"Down strong","Down")),"No-sig")
table(s2.phen$group1)
# Down Down strong      No-sig          Up   Up strong 
# 616         160       10030         634         717  
gene.s2 <- fread("/home/yuss/flyG4/result/Daniel.Robert.Genetics.RNAseq/002.1.gene.s2.bed") %>% as.data.frame()
s2.phen$num <- gene.s2[match(s2.phen$geneid,gene.s2$V4),7]
s2.phen$chr <- gene.s2[match(s2.phen$geneid,gene.s2$V4),1]
s2.phen <- s2.phen[s2.phen$chr=="X"|s2.phen$chr=="2L"|s2.phen$chr=="2R"|s2.phen$chr=="3L"|s2.phen$chr=="3R",]
s2.phen$chr1 <- ifelse(s2.phen$chr=="X","X","A")
s2.phen$start <- gene.s2[match(s2.phen$geneid,gene.s2$V4),2]
s2.phen$end <- gene.s2[match(s2.phen$geneid,gene.s2$V4),3]
s2.phen$length <- s2.phen$end-s2.phen$start
s2.phen$density <- s2.phen$num/s2.phen$length

s2.phen$group <- factor(s2.phen$group,levels = c("No-sig","Down","Up"))
my_comparisons = list(c("No-sig","Down"),c("No-sig","Up"),c("Down","Up"))

ggplot(data = s2.phen,aes(x=group,y=density*1000,fill=group)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(-0.18,-0.15,-0.175),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(name = '',
                    values = c('Up'='#D6604D','No-sig'='grey','Down'='#74ADD1'), #手动设置颜色时调整颜色的因子顺序
  )+
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Number of eG4 (per kb)") +
  coord_cartesian(ylim = c(0,0.3))

tapply(s2.phen$density, s2.phen$group, mean)

ggplot(data = s2.phen,aes(x=group,y=density*1000,fill=group)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  facet_grid(~chr1) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(-0.16,-0.13,-0.155),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(name = '',
                    values = c('Up'='#D6604D','No-sig'='grey','Down'='#74ADD1'), #手动设置颜色时调整颜色的因子顺序
  )+
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Number of eG4 (per kb)") +
  coord_cartesian(ylim = c(0,0.3))

s2.phen$group1 <- factor(s2.phen$group1,levels = c("No-sig","Down","Down strong","Up","Up strong"))
my_comparisons = list(c("No-sig","Down"),c("No-sig","Down strong"),c("No-sig","Up"),c("No-sig","Up strong"))

ggplot(data = s2.phen,aes(x=group1,y=density*1000,fill=group1)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  # facet_grid(~chr1) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(-0.1,-0.05,0,0.05),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(name = '',
                    values = c('Up'='#D6604D','No-sig'='grey','Down'='#74ADD1','Up strong'='#ae3b28','Down strong'='#3e8abb'), #手动设置颜色时调整颜色的因子顺序
  )+
  cowplot::theme_half_open()+
  theme(axis.text.x = element_text(size = 14,angle = 90, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Number of eG4 (per kb)") +
  coord_cartesian(ylim = c(0,0.5))

ggplot(data = s2.phen,aes(x=group1,y=density*1000,fill=group1)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  facet_grid(~chr1) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(-0.1,-0.05,0,0.05),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(name = '',
                    values = c('Up'='#D6604D','No-sig'='grey','Down'='#74ADD1','Up strong'='#ae3b28','Down strong'='#3e8abb'), #手动设置颜色时调整颜色的因子顺序
  )+
  cowplot::theme_half_open()+
  theme(axis.text.x = element_text(size = 14,angle = 90, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Number of eG4 (per kb)") +
  coord_cartesian(ylim = c(0,0.5))


s2.phen$numtype <- ifelse(s2.phen$num > 0,"have","no")
a <- as.data.frame(table(s2.phen$group, s2.phen$numtype))
wide_a <- spread(a, Var2, Freq)
wide_a$sum <- wide_a$have + wide_a$no
wide_a$have.ratio <- wide_a$have/wide_a$sum
wide_a$no.ratio <- wide_a$no/wide_a$sum


ggplot(wide_a, aes(x=Var1,y=have.ratio*100,fill=Var1)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(0.7),width = 0.7,color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值
  coord_cartesian(ylim = c(0,41)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_npg()+ ##更改颜色
  # geom_text(aes(label=have.ratio*100,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题(左下边框，没有网格线)
  ylab("(%) Gene with eG4") + ##通过bquote函数给图标签添加上下标
  geom_signif(y_position=35, xmin=1, xmax=2,
              annotation=c("6.0e-14"),tip_length=0)+ 
  geom_signif(y_position=37, xmin=2, xmax=3,
              annotation=c("0.0076"),tip_length=0)+ 
  geom_signif(y_position=39, xmin=1, xmax=3,
              annotation=c("5.5e-08"),tip_length=0)+ 
  theme(plot.title = element_text(size = 14),
        axis.title.y = element_text(size = 16), ##y坐标轴标题字体大小
        axis.title.x = element_blank(), ##删除x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.title = element_blank()) +
  labs(title="S2 cell")

data <- matrix((c(2118,7807,257,507)),nrow = 2)
fisher.test(data) # p-value = 5.956e-14
data <- matrix((c(2118,7807,377,967)),nrow = 2)
fisher.test(data) # p-value =5.491e-08
data <- matrix((c(377,967,257,507)),nrow = 2)
fisher.test(data) # p-value = 0.007644
#group1
s2.phen$numtype <- ifelse(s2.phen$num > 0,"have","no")
a <- as.data.frame(table(s2.phen$group1, s2.phen$numtype))
wide_a <- spread(a, Var2, Freq)
wide_a$sum <- wide_a$have + wide_a$no
wide_a$have.ratio <- wide_a$have/wide_a$sum
wide_a$no.ratio <- wide_a$no/wide_a$sum

ggplot(wide_a, aes(x=Var1,y=have.ratio*100,fill=Var1)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(0.7),width = 0.7,color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值
  coord_cartesian(ylim = c(0,50)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_npg()+ ##更改颜色
  # geom_text(aes(label=have.ratio*100,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题(左下边框，没有网格线)
  ylab("(%) Gene with eG4") + ##通过bquote函数给图标签添加上下标
  geom_signif(y_position=48, xmin=3, xmax=5,
              annotation=c("***"),tip_length=0)+ 
  geom_signif(y_position=44, xmin=3, xmax=4,
              annotation=c("**"),tip_length=0)+ 
  theme(plot.title = element_text(size = 14),
        axis.title.y = element_text(size = 16), ##y坐标轴标题字体大小
        axis.title.x = element_blank(), ##删除x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.title = element_blank()) +
  labs(title="S2 cell")

#chrX,chA
s2.phen$numtype <- ifelse(s2.phen$num > 0,"have","no")
a <- as.data.frame(table(s2.phen$group,s2.phen$chr1, s2.phen$numtype))
wide_a <- spread(a, Var3, Freq)
wide_a$sum <- wide_a$have + wide_a$no
wide_a$have.ratio <- wide_a$have/wide_a$sum
wide_a$no.ratio <- wide_a$no/wide_a$sum

ggplot(wide_a, aes(x=Var1,y=have.ratio*100,fill=Var2)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(0.7),width = 0.7,color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值
  coord_cartesian(ylim = c(0,41)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_npg()+ ##更改颜色
  # geom_text(aes(label=have.ratio*100,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题(左下边框，没有网格线)
  ylab("(%) Gene with eG4") + ##通过bquote函数给图标签添加上下标
  geom_signif(y_position=35, xmin=1.2, xmax=2.2,
              annotation=c("0.017"),tip_length=0)+ 
  geom_signif(y_position=37, xmin=2.2, xmax=3.2,
              annotation=c("0.89"),tip_length=0)+ 
  geom_signif(y_position=39, xmin=1.2, xmax=3.2,
              annotation=c("1.2e-04"),tip_length=0)+ 
  theme(plot.title = element_text(size = 14),
        axis.title.y = element_text(size = 16), ##y坐标轴标题字体大小
        axis.title.x = element_blank(), ##删除x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.title = element_blank()) +
  labs(title="S2 cell")


data <- matrix((c(266,1337,23,62)),nrow = 2)
fisher.test(data) # 0.01745

data <- matrix((c(266,1337,83,234)),nrow = 2)
fisher.test(data) # p-value = 0.0001196

data <- matrix((c(23,62,83,234)),nrow = 2)
fisher.test(data) # p-value = 0.0001196

#chrX,chA group1
s2.phen$numtype <- ifelse(s2.phen$num > 0,"have","no")
a <- as.data.frame(table(s2.phen$group1,s2.phen$chr1, s2.phen$numtype))
wide_a <- spread(a, Var3, Freq)
wide_a$sum <- wide_a$have + wide_a$no
wide_a$have.ratio <- wide_a$have/wide_a$sum
wide_a$no.ratio <- wide_a$no/wide_a$sum

ggplot(wide_a, aes(x=Var1,y=have.ratio*100,fill=Var2)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(0.7),width = 0.7,color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值
  coord_cartesian(ylim = c(0,50)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_npg()+ ##更改颜色
  # geom_text(aes(label=have.ratio*100,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题(左下边框，没有网格线)
  ylab("(%) Gene with eG4") + ##通过bquote函数给图标签添加上下标
  # geom_signif(y_position=28, xmin=2.2, xmax=3.2,
  #             annotation=c("****"),tip_length=0)+ #p-value = 0.4004
  theme(plot.title = element_text(size = 14),
        axis.title.y = element_text(size = 16), ##y坐标轴标题字体大小
        axis.title.x = element_blank(), ##删除x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.title = element_blank()) +
  labs(title="S2 cell")

#-----------filter-------------------------
# 不能把0过滤掉，如果把0剔除，那算出来的密度也不准确

#### 差异基因分成四类看启动子S2 G4的数量 ####
rm(list = ls());gc();rm(list = ls())#清空
Num = "010.3."
promoter.s2G4 <- fread("/home/yuss/flyG4/result/LucyCherbas.GR.2010.RNAseq/007.4.promoter2000.s2.bed") %>% as.data.frame()
s2.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.DE.s2.Phen.txt") %>% as.data.frame()
s2.Phen$group1 <- ifelse(s2.Phen$pvalue<0.05&abs(s2.Phen$log2FoldChange)>=0.5,ifelse(s2.Phen$log2FoldChange>0.5,ifelse(s2.Phen$pvalue<0.01&s2.Phen$log2FoldChange>1,"Up strong","Up"),ifelse(s2.Phen$pvalue<0.01&s2.Phen$log2FoldChange<=-1,"Down strong","Down")),"No-sig")
table(s2.Phen$group1)
s2.Phen$promoter.num <- promoter.s2G4[match(s2.Phen$geneid,promoter.s2G4$V4),6]
s2.Phen$chr <- promoter.s2G4[match(s2.Phen$geneid,promoter.s2G4$V4),1]
s2.Phen <- s2.Phen[s2.Phen$chr=="X"|s2.Phen$chr=="2L"|s2.Phen$chr=="2R"|s2.Phen$chr=="3L"|s2.Phen$chr=="3R",]
s2.Phen$chr1 <- ifelse(s2.Phen$chr=="X","X","A")
s2.Phen$group <- factor(s2.Phen$group,levels = c("No-sig","Down","Up"))
my_comparisons = list(c("No-sig","Down"),c("No-sig","Up"))
ggplot(data = s2.Phen,aes(x=group,y=promoter.num,fill=group)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(3,3.5,4,4.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  # scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Number of Promoter eG4") +
  coord_cartesian(ylim = c(0,10)) +
  scale_y_continuous(expand = c(0,0)) 
tapply(s2.Phen$promoter.num, s2.Phen$group, mean)

ggplot(data = s2.Phen,aes(x=group,y=promoter.num,fill=group)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  facet_grid(~chr1)+
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(3,3.5,4,4.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  # scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Number of Promoter eG4") +
  coord_cartesian(ylim = c(0,10)) +
  scale_y_continuous(expand = c(0,0)) 
tapply(s2.Phen$promoter.num, s2.Phen$group, mean)
aggregate(s2.Phen$promoter.num, by = list(s2.Phen$group, s2.Phen$chr1), FUN = median)

s2.Phen$group1 <- factor(s2.Phen$group1,levels = c("No-sig","Down","Down strong","Up","Up strong"))
my_comparisons = list(c("No-sig","Down"),c("No-sig","Down strong"),c("No-sig","Up"),c("No-sig","Up strong"))
ggplot(data = s2.Phen,aes(x=group1,y=promoter.num,fill=group1)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(3,3.5,4,4.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  # scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Number of Promoter eG4") +
  coord_cartesian(ylim = c(0,5))# +
scale_y_continuous(expand = c(0,0)) 



ggplot(data = s2.Phen,aes(x=group1,y=promoter.num,fill=group1)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  facet_grid(~chr1) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(3,3.5,4,4.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  # scale_fill_manual(values = c("#4d4d4d","#FDDBC7","#92c5de")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Number of Promoter eG4") +
  coord_cartesian(ylim = c(0,5))# +
scale_y_continuous(expand = c(0,0)) 

##把0过滤
s2.Phen.filter <- s2.Phen[s2.Phen$promoter.num != 0, ]
s2.Phen.filter$group1 <- factor(s2.Phen.filter$group1,levels = c("No-sig","Down","Down strong","Up","Up strong"))
my_comparisons = list(c("No-sig","Down"),c("No-sig","Down strong"),c("No-sig","Up"),c("No-sig","Up strong"))
ggplot(data = s2.Phen.filter,aes(x=group1,y=promoter.num,fill=group1)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(3,3.5,4,4.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(name = '',
                    values = c('Up'='#D6604D','No-sig'='grey','Down'='#74ADD1','Up strong'='#ae3b28','Down strong'='#3e8abb') )+#手动设置颜色时调整颜色的因子顺序
  cowplot::theme_half_open()+
  theme(axis.text.x = element_text(size = 14,angle = 90, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Number of Promoter eG4") +
  coord_cartesian(ylim = c(0,5)) +
  scale_y_continuous(expand = c(0,0)) 

tapply(s2.Phen.filter$promoter.num, s2.Phen.filter$group1, mean)


ggplot(data = s2.Phen.filter,aes(x=group1,y=promoter.num,fill=group1)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  facet_grid(~chr1) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(3,3.5,4,4.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(name = '',
                    values = c('Up'='#D6604D','No-sig'='grey','Down'='#74ADD1','Up strong'='#ae3b28','Down strong'='#3e8abb') )+#手动设置颜色时调整颜色的因子顺序
  cowplot::theme_half_open()+
  theme(axis.text.x = element_text(size = 14,angle = 90, hjust = 0.5),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Number of Promoter eG4") +
  coord_cartesian(ylim = c(0,5)) +
  scale_y_continuous(expand = c(0,0)) 

tapply(s2.Phen.filter$promoter.num, s2.Phen.filter$group1, mean)


s2.Phen.filter$group <- factor(s2.Phen.filter$group,levels = c("No-sig","Down","Up"))
my_comparisons = list(c("No-sig","Down"),c("No-sig","Up"))
ggplot(data = s2.Phen.filter,aes(x=group,y=promoter.num,fill=group)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(3,3.5,4,4.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(name = '',
                    values = c('Up'='#D6604D','No-sig'='grey','Down'='#74ADD1'), #手动设置颜色时调整颜色的因子顺序
  )+
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Number of Promoter eG4") +
  coord_cartesian(ylim = c(0,5)) +
  scale_y_continuous(expand = c(0,0))

tapply(s2.Phen.filter$promoter.num, s2.Phen.filter$group, mean)

ggplot(data = s2.Phen.filter,aes(x=group,y=promoter.num,fill=group)) +
  geom_boxplot(notch = F,outlier.colour = "white") +
  facet_grid(~chr1) +
  stat_compare_means(comparisons = my_comparisons,
                     label.y = c(3,3.5,4,4.5),
                     aes(label = paste0("p = ", ..p.format..)),tip.length = 0,size=4)+
  scale_fill_manual(name = '',
                    values = c('Up'='#D6604D','No-sig'='grey','Down'='#74ADD1'), #手动设置颜色时调整颜色的因子顺序
  )+
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("Number of eG4 in Promoter") +
  coord_cartesian(ylim = c(0,5)) +
  scale_y_continuous(expand = c(0,0))

s2.Phen$promter.numtype <- ifelse(s2.Phen$promoter.num > 0,"have","no")
a <- as.data.frame(table(s2.Phen$group, s2.Phen$promter.numtype))
wide_a <- spread(a, Var2, Freq)
wide_a$sum <- wide_a$have + wide_a$no
wide_a$have.ratio <- wide_a$have/wide_a$sum
wide_a$no.ratio <- wide_a$no/wide_a$sum


ggplot(wide_a, aes(x=Var1,y=have.ratio*100,fill=Var1)) + ##fill是图形的填充色
  geom_bar(stat = 'identity',position = position_dodge(0.7),width = 0.7,color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值
  coord_cartesian(ylim = c(0,40)) + ##坐标轴范围
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_npg()+ ##更改颜色
  # geom_text(aes(label=have.ratio*100,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
  cowplot::theme_half_open() + ##主题(左下边框，没有网格线)
  ylab("(%) Gene with Promoter eG4") + ##通过bquote函数给图标签添加上下标
  geom_signif(y_position=25, xmin=1, xmax=2,
              annotation=c("0.964"),tip_length=0)+
  geom_signif(y_position=35, xmin=1, xmax=3,
              annotation=c("5.1e-06"),tip_length=0)+
  geom_signif(y_position=31, xmin=2, xmax=3,
              annotation=c("0.004"),tip_length=0)+
  theme(plot.title = element_text(size = 14),
        axis.title.y = element_text(size = 16), ##y坐标轴标题字体大小
        axis.title.x = element_blank(), ##删除x坐标轴标题
        axis.text = element_text(size=14), ##轴文本字体大小
        legend.title = element_blank()) +
  labs(title="S2 cell")

data <- matrix((c(2268,7657,175,589)),nrow = 2)
fisher.test(data) # p-value = 0.9644
data <- matrix((c(2268,7657,384,960)),nrow = 2)
fisher.test(data) # p-value = 5.103e-06
data <- matrix((c(384,960,175,589)),nrow = 2)
fisher.test(data) # p-value = 0.004752
#chrX,所有常染色体

  #chrX,chA
  s2.Phen$promter.numtype <- ifelse(s2.Phen$promoter.num > 0,"have","no")
  a <- as.data.frame(table(s2.Phen$group,s2.Phen$chr1, s2.Phen$promter.numtype))
  wide_a <- spread(a, Var3, Freq)
  wide_a$sum <- wide_a$have + wide_a$no
  wide_a$have.ratio <- wide_a$have/wide_a$sum
  wide_a$no.ratio <- wide_a$no/wide_a$sum

ggplot(wide_a, aes(x=Var1,y=have.ratio*100,fill=Var2)) + ##fill是图形的填充色
    geom_bar(stat = 'identity',position = position_dodge(0.7),width = 0.7,color = 'white') + ##stat：设置统计方法,identity表示条形的高度是变量的值
    coord_cartesian(ylim = c(0,40)) + ##坐标轴范围
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_npg()+ ##更改颜色
    # geom_text(aes(label=have.ratio*100,vjust = -0.5),color="black", size=4) + ##柱形图上加数值标签
    cowplot::theme_half_open() + ##主题(左下边框，没有网格线)
    ylab("(%) Gene with Promoter eG4") + ##通过bquote函数给图标签添加上下标
  geom_signif(y_position=35, xmin=1.2, xmax=3.2,
              annotation=c("0.004"),tip_length=0)+ #p-value = 0.4004
    theme(plot.title = element_text(size = 14),
          axis.title.y = element_text(size = 16), ##y坐标轴标题字体大小
          axis.title.x = element_blank(), ##删除x坐标轴标题
          axis.text = element_text(size=14), ##轴文本字体大小
          legend.title = element_blank()) +
    labs(title="S2 cell")

data <- matrix((c(275,1328,17,68)),nrow = 2)
fisher.test(data) # p-value = 0.465
  
data <- matrix((c(275,1328,77,240)),nrow = 2)
fisher.test(data) # p-value =0.004126
data <- matrix((c(17,68,77,240)),nrow = 2)
fisher.test(data) # p-value =0.4717
