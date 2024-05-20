rm(list = ls());gc();rm(list = ls())
Num = "001.11."
#*PhyloP score------------------------------------------------------------------------------------
phyloop1G <- fread("/home/yuss/flyG4/result/PQS/001.11.phyloP_loop1.includeG.bed") %>% data.frame()
phyloop2G <- fread("/home/yuss/flyG4/result/PQS/001.11.phyloP_loop2.includeG.bed") %>% data.frame()
phyloop3G <- fread("/home/yuss/flyG4/result/PQS/001.11.phyloP_loop3.includeG.bed") %>% data.frame()
phyloop1NG <- fread("/home/yuss/flyG4/result/PQS/001.11.phyloP_loop1.NincludeG.bed") %>% data.frame()
phyloop2NG <- fread("/home/yuss/flyG4/result/PQS/001.11.phyloP_loop2.NincludeG.bed") %>% data.frame()
phyloop3NG <- fread("/home/yuss/flyG4/result/PQS/001.11.phyloP_loop3.NincludeG.bed") %>% data.frame()

merge.all <- fread("/home/yuss/flyG4/result/PQS/001.2.merge.all.bed") %>% data.frame()
merge.all$phyloop1G <- phyloop1G[match(merge.all$id,phyloop1G$V1),6]
merge.all$phyloop2G <- phyloop2G[match(merge.all$id,phyloop2G$V1),6]
merge.all$phyloop3G <- phyloop3G[match(merge.all$id,phyloop3G$V1),6]
df.phyloop1G <- na.omit(merge.all[,c("type","phyloop1G")])
df.phyloop1G$class <- "include_G"
df.phyloop2G <- na.omit(merge.all[,c("type","phyloop2G")])
df.phyloop2G$class <- "include_G"
df.phyloop3G <- na.omit(merge.all[,c("type","phyloop3G")])
df.phyloop3G$class <- "include_G"
merge.all$phyloop1NG <- phyloop1NG[match(merge.all$id,phyloop1NG$V1),6]
merge.all$phyloop2NG <- phyloop2NG[match(merge.all$id,phyloop2NG$V1),6]
merge.all$phyloop3NG <- phyloop3NG[match(merge.all$id,phyloop3NG$V1),6]
df.phyloop1NG <- na.omit(merge.all[,c("type","phyloop1NG")])
df.phyloop1NG$class <- "without_G"
df.phyloop2NG <- na.omit(merge.all[,c("type","phyloop2NG")])
df.phyloop2NG$class <- "without_G"
df.phyloop3NG <- na.omit(merge.all[,c("type","phyloop3NG")])
df.phyloop3NG$class <- "without_G"

##修改列名，列名不一样不能合并
colnames(df.phyloop1G)[2] <- "phy"
colnames(df.phyloop2G)[2] <- "phy"
colnames(df.phyloop3G)[2] <- "phy"
colnames(df.phyloop1NG)[2] <- "phy"
colnames(df.phyloop2NG)[2] <- "phy"
colnames(df.phyloop3NG)[2] <- "phy"
df <- rbind(df.phyloop1G,df.phyloop2G,df.phyloop3G,df.phyloop1NG,df.phyloop2NG,df.phyloop3NG)

##plot
df$type <- factor(df$type,levels = c("non_eG4","kc_specific","s2_specific","overlap","merge"))
library(ggpubr) 
ggplot(data = df,aes(x=type,y=phy,fill=class)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  scale_fill_manual(values = c("#fdb863","#74ADD1")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank()) +
  ylab("PhyloP score") +
  stat_compare_means(aes(group = class),method = "t.test",label = "p.signif",label.y = 3) +
  coord_cartesian(ylim = c(-2,3)) 
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"G4IncludeOrWithoutGLoop.PhyloP.pdf"),
       device = "pdf",width = 7,height = 5)

#*PhastCons score------------------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())
Num = "001.11."
##读取数据
phaloop1G <- fread("/home/yuss/flyG4/result/PQS/001.11.phastCons_loop1.includeG.bed") %>% data.frame()
phaloop2G <- fread("/home/yuss/flyG4/result/PQS/001.11.phastCons_loop2.includeG.bed") %>% data.frame()
phaloop3G <- fread("/home/yuss/flyG4/result/PQS/001.11.phastCons_loop3.includeG.bed") %>% data.frame()
phaloop1NG <- fread("/home/yuss/flyG4/result/PQS/001.11.phastCons_loop1.NincludeG.bed") %>% data.frame()
phaloop2NG <- fread("/home/yuss/flyG4/result/PQS/001.11.phastCons_loop2.NincludeG.bed") %>% data.frame()
phaloop3NG <- fread("/home/yuss/flyG4/result/PQS/001.11.phastCons_loop3.NincludeG.bed") %>% data.frame()

merge.all <- fread("/home/yuss/flyG4/result/PQS/001.2.merge.all.bed") %>% data.frame()
merge.all$phaloop1G <- phaloop1G[match(merge.all$id,phaloop1G$V1),6]
merge.all$phaloop2G <- phaloop2G[match(merge.all$id,phaloop2G$V1),6]
merge.all$phaloop3G <- phaloop3G[match(merge.all$id,phaloop3G$V1),6]
df.phaloop1G <- na.omit(merge.all[,c("type","phaloop1G")])
df.phaloop1G$class <- "include_G"
df.phaloop2G <- na.omit(merge.all[,c("type","phaloop2G")])
df.phaloop2G$class <- "include_G"
df.phaloop3G <- na.omit(merge.all[,c("type","phaloop3G")])
df.phaloop3G$class <- "include_G"
merge.all$phaloop1NG <- phaloop1NG[match(merge.all$id,phaloop1NG$V1),6]
merge.all$phaloop2NG <- phaloop2NG[match(merge.all$id,phaloop2NG$V1),6]
merge.all$phaloop3NG <- phaloop3NG[match(merge.all$id,phaloop3NG$V1),6]
df.phaloop1NG <- na.omit(merge.all[,c("type","phaloop1NG")])
df.phaloop1NG$class <- "without_G"
df.phaloop2NG <- na.omit(merge.all[,c("type","phaloop2NG")])
df.phaloop2NG$class <- "without_G"
df.phaloop3NG <- na.omit(merge.all[,c("type","phaloop3NG")])
df.phaloop3NG$class <- "without_G"

##修改列名，列名不一样不能合并
colnames(df.phaloop1G)[2] <- "pha"
colnames(df.phaloop2G)[2] <- "pha"
colnames(df.phaloop3G)[2] <- "pha"
colnames(df.phaloop1NG)[2] <- "pha"
colnames(df.phaloop2NG)[2] <- "pha"
colnames(df.phaloop3NG)[2] <- "pha"

##合并dataframe
df <- rbind(df.phaloop1G,df.phaloop2G,df.phaloop3G,df.phaloop1NG,df.phaloop2NG,df.phaloop3NG)

##plot
df$type <- factor(df$type,levels = c("non_eG4","kc_specific","s2_specific","overlap","merge"))
library(ggpubr) 
ggplot(data = df,aes(x=type,y=pha,fill=class)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  scale_fill_manual(values = c("#fdb863","#74ADD1")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank()) +
  ylab("PhastCons score") +
  stat_compare_means(aes(group = class),method = "t.test",label = "p.signif",label.y = 1.15) +
  coord_cartesian(ylim = c(0,1.2)) 
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"G4IncludeOrWithoutGLoop.PhastCons.pdf"),
device = "pdf",width = 7,height = 5)


rm(list = ls());gc();rm(list = ls())
Num = "001.11."
#*剔除序列后的PhyloP score------------------------------------------------------------------------------------
phyloop1G <- fread("/home/yuss/flyG4/result/PQS/001.11.phyloP_loop1.includeG.delete0.bed") %>% data.frame()
phyloop2G <- fread("/home/yuss/flyG4/result/PQS/001.11.phyloP_loop2.includeG.delete0.bed") %>% data.frame()
phyloop3G <- fread("/home/yuss/flyG4/result/PQS/001.11.phyloP_loop3.includeG.delete0.bed") %>% data.frame()
phyloop1NG <- fread("/home/yuss/flyG4/result/PQS/001.11.phyloP_loop1.NincludeG.delete0.bed") %>% data.frame()
phyloop2NG <- fread("/home/yuss/flyG4/result/PQS/001.11.phyloP_loop2.NincludeG.delete0.bed") %>% data.frame()
phyloop3NG <- fread("/home/yuss/flyG4/result/PQS/001.11.phyloP_loop3.NincludeG.delete0.bed") %>% data.frame()

merge.all <- fread("/home/yuss/flyG4/result/PQS/001.2.merge.all.bed") %>% data.frame()
merge.all$phyloop1G <- phyloop1G[match(merge.all$id,phyloop1G$V1),6]
merge.all$phyloop2G <- phyloop2G[match(merge.all$id,phyloop2G$V1),6]
merge.all$phyloop3G <- phyloop3G[match(merge.all$id,phyloop3G$V1),6]
df.phyloop1G <- na.omit(merge.all[,c("type","phyloop1G")])
df.phyloop1G$class <- "include_G"
df.phyloop2G <- na.omit(merge.all[,c("type","phyloop2G")])
df.phyloop2G$class <- "include_G"
df.phyloop3G <- na.omit(merge.all[,c("type","phyloop3G")])
df.phyloop3G$class <- "include_G"
merge.all$phyloop1NG <- phyloop1NG[match(merge.all$id,phyloop1NG$V1),6]
merge.all$phyloop2NG <- phyloop2NG[match(merge.all$id,phyloop2NG$V1),6]
merge.all$phyloop3NG <- phyloop3NG[match(merge.all$id,phyloop3NG$V1),6]
df.phyloop1NG <- na.omit(merge.all[,c("type","phyloop1NG")])
df.phyloop1NG$class <- "without_G"
df.phyloop2NG <- na.omit(merge.all[,c("type","phyloop2NG")])
df.phyloop2NG$class <- "without_G"
df.phyloop3NG <- na.omit(merge.all[,c("type","phyloop3NG")])
df.phyloop3NG$class <- "without_G"

##修改列名，列名不一样不能合并
colnames(df.phyloop1G)[2] <- "phy"
colnames(df.phyloop2G)[2] <- "phy"
colnames(df.phyloop3G)[2] <- "phy"
colnames(df.phyloop1NG)[2] <- "phy"
colnames(df.phyloop2NG)[2] <- "phy"
colnames(df.phyloop3NG)[2] <- "phy"
df <- rbind(df.phyloop1G,df.phyloop2G,df.phyloop3G,df.phyloop1NG,df.phyloop2NG,df.phyloop3NG)

##plot
df$type <- factor(df$type,levels = c("non_eG4","kc_specific","s2_specific","overlap","merge"))
library(ggpubr) 
ggplot(data = df,aes(x=type,y=phy,fill=class)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  scale_fill_manual(values = c("#fdb863","#74ADD1")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank()) +
  ylab("PhyloP score") +
  stat_compare_means(aes(group = class),method = "t.test",label = "p.signif",label.y = 3) +
  coord_cartesian(ylim = c(-2,3)) 
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"G4IncludeOrWithoutGLoop.PhyloP.pdf"),
       device = "pdf",width = 7,height = 5)


#*剔除序列后的PhastCons score------------------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())
Num = "001.11."
##读取数据
phaloop1G <- fread("/home/yuss/flyG4/result/PQS/001.11.phastCons_loop1.includeG.bed") %>% data.frame()
phaloop2G <- fread("/home/yuss/flyG4/result/PQS/001.11.phastCons_loop2.includeG.bed") %>% data.frame()
phaloop3G <- fread("/home/yuss/flyG4/result/PQS/001.11.phastCons_loop3.includeG.bed") %>% data.frame()
phaloop1NG <- fread("/home/yuss/flyG4/result/PQS/001.11.phastCons_loop1.NincludeG.delete-nan.bed") %>% data.frame()
phaloop2NG <- fread("/home/yuss/flyG4/result/PQS/001.11.phastCons_loop2.NincludeG.delete-nan.bed") %>% data.frame()
phaloop3NG <- fread("/home/yuss/flyG4/result/PQS/001.11.phastCons_loop3.NincludeG.delete-nan.bed") %>% data.frame()

merge.all <- fread("/home/yuss/flyG4/result/PQS/001.2.merge.all.bed") %>% data.frame()
merge.all$phaloop1G <- phaloop1G[match(merge.all$id,phaloop1G$V1),6]
merge.all$phaloop2G <- phaloop2G[match(merge.all$id,phaloop2G$V1),6]
merge.all$phaloop3G <- phaloop3G[match(merge.all$id,phaloop3G$V1),6]
df.phaloop1G <- na.omit(merge.all[,c("type","phaloop1G")])
df.phaloop1G$class <- "include_G"
df.phaloop2G <- na.omit(merge.all[,c("type","phaloop2G")])
df.phaloop2G$class <- "include_G"
df.phaloop3G <- na.omit(merge.all[,c("type","phaloop3G")])
df.phaloop3G$class <- "include_G"
merge.all$phaloop1NG <- phaloop1NG[match(merge.all$id,phaloop1NG$V1),6]
merge.all$phaloop2NG <- phaloop2NG[match(merge.all$id,phaloop2NG$V1),6]
merge.all$phaloop3NG <- phaloop3NG[match(merge.all$id,phaloop3NG$V1),6]
df.phaloop1NG <- na.omit(merge.all[,c("type","phaloop1NG")])
df.phaloop1NG$class <- "without_G"
df.phaloop2NG <- na.omit(merge.all[,c("type","phaloop2NG")])
df.phaloop2NG$class <- "without_G"
df.phaloop3NG <- na.omit(merge.all[,c("type","phaloop3NG")])
df.phaloop3NG$class <- "without_G"

##修改列名，列名不一样不能合并
colnames(df.phaloop1G)[2] <- "pha"
colnames(df.phaloop2G)[2] <- "pha"
colnames(df.phaloop3G)[2] <- "pha"
colnames(df.phaloop1NG)[2] <- "pha"
colnames(df.phaloop2NG)[2] <- "pha"
colnames(df.phaloop3NG)[2] <- "pha"

##合并dataframe
df <- rbind(df.phaloop1G,df.phaloop2G,df.phaloop3G,df.phaloop1NG,df.phaloop2NG,df.phaloop3NG)

##plot
df$type <- factor(df$type,levels = c("non_eG4","kc_specific","s2_specific","overlap","merge"))
library(ggpubr) 
ggplot(data = df,aes(x=type,y=pha,fill=class)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  scale_fill_manual(values = c("#fdb863","#74ADD1")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank()) +
  ylab("PhastCons score") +
  stat_compare_means(aes(group = class),method = "t.test",label = "p.signif",label.y = 1.15) +
  coord_cartesian(ylim = c(0,1.2)) 
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"G4IncludeOrWithoutGLoop.PhastCons.pdf"),
       device = "pdf",width = 7,height = 5)