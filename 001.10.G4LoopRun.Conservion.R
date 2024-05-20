rm(list = ls());gc();rm(list = ls())#清空
Num = "001.10."

pqs <- fread("/home/yuss/flyG4/result/PQS/001.10.dmel.pqs.chr.bed") %>% as.data.frame()
colnames(pqs) <- c("chr","start","end","width","id","rl1","rl2","rl3","ll1","ll2","ll3")
pqs$start <- pqs$start - 1
run1 <- pqs[,c(1,2,5)]
run1$end <- pqs$start + pqs$rl1
run1 <- run1[,c(1,2,4,3)]

loop1 <- pqs[,c(1,5)]
loop1$start <- run1$end
loop1$end <- loop1$start + pqs$ll1
loop1 <- loop1[,c(1,3,4,2)]

run2 <- pqs[,c(1,5)]
run2$start <- loop1$end
run2$end <- run2$start + pqs$rl2
run2 <- run2[,c(1,3,4,2)]

loop2 <- pqs[,c(1,5)]
loop2$start <- run2$end
loop2$end <- loop2$start + pqs$ll2
loop2 <- loop2[,c(1,3,4,2)]

run3 <- pqs[,c(1,5)]
run3$start <- loop2$end
run3$end <- run3$start + pqs$rl3
run3 <- run3[,c(1,3,4,2)]

loop3 <- pqs[,c(1,5)]
loop3$start <- run3$end
loop3$end <- loop3$start + pqs$ll3
loop3 <- loop3[,c(1,3,4,2)]

run4 <- pqs[,c(1,3,5)]
run4$start <- loop3$end
run4 <- run4[,c(1,4,2,3)]
write.table(run1,file = paste0("/home/yuss/flyG4/result/PQS/",Num,"pqs_run1.bed"),sep = '\t',
            col.names = F,row.names = F,quote = F)
write.table(run2,file = paste0("/home/yuss/flyG4/result/PQS/",Num,"pqs_run2.bed"),sep = '\t',
            col.names = F,row.names = F,quote = F)
write.table(run3,file = paste0("/home/yuss/flyG4/result/PQS/",Num,"pqs_run3.bed"),sep = '\t',
            col.names = F,row.names = F,quote = F)
write.table(run4,file = paste0("/home/yuss/flyG4/result/PQS/",Num,"pqs_run4.bed"),sep = '\t',
            col.names = F,row.names = F,quote = F)
write.table(loop1,file = paste0("/home/yuss/flyG4/result/PQS/",Num,"pqs_loop1.bed"),sep = '\t',
            col.names = F,row.names = F,quote = F)
write.table(loop2,file = paste0("/home/yuss/flyG4/result/PQS/",Num,"pqs_loop2.bed"),sep = '\t',
            col.names = F,row.names = F,quote = F)
write.table(loop3,file = paste0("/home/yuss/flyG4/result/PQS/",Num,"pqs_loop3.bed"),sep = '\t',
            col.names = F,row.names = F,quote = F)

#*phastCons--------------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())#清空
Num = "001.10."
setwd("/home/yuss/flyG4/result/PQS")
filename <- list.files("/home/yuss/flyG4/result/PQS",pattern = "^001.10.phastCons")
for (i in 1:length(filename)){
  var_name <- gsub('.bed', '',filename[i])
  var_name <- gsub('001.10.phastCons_', 'pha_',var_name)
  assign(var_name, read.table(filename[i], sep = '\t', header = F)) ##assign()函数将一个读取的数据框对象分配给先前定义的变量名 var_name
}

merge.all <- fread("/home/yuss/flyG4/result/PQS/001.2.merge.all.bed") %>% data.frame()
merge.all$pha_loop1 <- pha_loop1[match(merge.all$id,pha_loop1$V1),6]
merge.all$pha_loop2 <- pha_loop2[match(merge.all$id,pha_loop2$V1),6]
merge.all$pha_loop3 <- pha_loop3[match(merge.all$id,pha_loop3$V1),6]
merge.all$pha_run1 <- pha_run1[match(merge.all$id,pha_run1$V1),6]
merge.all$pha_run2 <- pha_run2[match(merge.all$id,pha_run2$V1),6]
merge.all$pha_run3 <- pha_run3[match(merge.all$id,pha_run3$V1),6]
merge.all$pha_run4 <- pha_run4[match(merge.all$id,pha_run4$V1),6]
df_pha_loop1 <- merge.all[,c("type","pha_loop1")]
df_pha_loop1$class <- "loop"
colnames(df_pha_loop1)[2] <- "pha"
df_pha_loop2 <- merge.all[,c("type","pha_loop2")]
df_pha_loop2$class <- "loop"
colnames(df_pha_loop2)[2] <- "pha"
df_pha_loop3 <- merge.all[,c("type","pha_loop3")]
df_pha_loop3$class <- "loop"
colnames(df_pha_loop3)[2] <- "pha"
df_pha_run1 <- merge.all[,c("type","pha_run1")]
df_pha_run1$class <- "run"
colnames(df_pha_run1)[2] <- "pha"
df_pha_run2 <- merge.all[,c("type","pha_run2")]
df_pha_run2$class <- "run"
colnames(df_pha_run2)[2] <- "pha"
df_pha_run3 <- merge.all[,c("type","pha_run3")]
df_pha_run3$class <- "run"
colnames(df_pha_run3)[2] <- "pha"
df_pha_run4 <- merge.all[,c("type","pha_run4")]
df_pha_run4$class <- "run"
colnames(df_pha_run4)[2] <- "pha"

data_frames_to_merge <- list(
  df_pha_loop1,
  df_pha_loop2,
  df_pha_loop3,
  df_pha_run1,
  df_pha_run2,
  df_pha_run3,
  df_pha_run4
)

# 使用 Reduce 函数按列合并数据框
df <- Reduce(rbind, data_frames_to_merge)
df$type <- factor(df$type,levels =c("non_eG4","kc_specific","s2_specific","overlap","merge"))

# 查看表格中是否有缺失值
if (anyNA(df)) {
  cat("表格中存在缺失值。\n")
} else {
  cat("表格中没有缺失值。\n")
}

ggplot(data = df,aes(x=type,y=pha,fill=class)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  scale_fill_manual(values = c("#74ADD1","#D6604D")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank()) +
  ylab("PhastCons score") +
  stat_compare_means(aes(group = class),method = "t.test",label = "p.signif",label.y = 1.15) +
  coord_cartesian(ylim = c(0,1.2)) 
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"G4LoopRun.PhastCons.pdf"),
         device = "pdf",width = 7,height = 5)


#*phyloP--------------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())#清空
Num = "001.10."

setwd("/home/yuss/flyG4/result/PQS")
filename <- list.files("/home/yuss/flyG4/result/PQS",pattern = "^001.10.phyloP")
for (i in 1:length(filename)){
  var_name <- gsub('.bed', '',filename[i])
  var_name <- gsub('001.10.phyloP_', 'phy_',var_name)
  assign(var_name, read.table(filename[i], sep = '\t', header = F)) ##assign()函数将一个读取的数据框对象分配给先前定义的变量名 var_name
}

merge.all <- fread("/home/yuss/flyG4/result/PQS/001.2.merge.all.bed") %>% data.frame()

# 假设您有一个列表，其中包含了 pha_loop1、pha_loop2、pha_loop3、pha_run1 等数据框
# 这里假设这些数据框按照一定顺序存储在一个列表中
data.frames <- list(phy_loop1, phy_loop2, phy_loop3, phy_run1, phy_run2, phy_run3, phy_run4)

# 定义要创建的新列名
new.colnames <- c("phy_loop1", "phy_loop2", "phy_loop3", "phy_run1", "phy_run2", "phy_run3", "phy_run4")

# 循环操作
for (i in seq_along(data.frames)) {
  colname <- new.colnames[i]
  merge.all[[colname]] <- data.frames[[i]][match(merge.all$id, data.frames[[i]]$V1), 6]
}

df_phy_loop1 <- merge.all[,c("type","phy_loop1")]
df_phy_loop1$class <- "loop"
colnames(df_phy_loop1)[2] <- "phy"
df_phy_loop2 <- merge.all[,c("type","phy_loop2")]
df_phy_loop2$class <- "loop"
colnames(df_phy_loop2)[2] <- "phy"
df_phy_loop3 <- merge.all[,c("type","phy_loop3")]
df_phy_loop3$class <- "loop"
colnames(df_phy_loop3)[2] <- "phy"
df_phy_run1 <- merge.all[,c("type","phy_run1")]
df_phy_run1$class <- "run"
colnames(df_phy_run1)[2] <- "phy"
df_phy_run2 <- merge.all[,c("type","phy_run2")]
df_phy_run2$class <- "run"
colnames(df_phy_run2)[2] <- "phy"
df_phy_run3 <- merge.all[,c("type","phy_run3")]
df_phy_run3$class <- "run"
colnames(df_phy_run3)[2] <- "phy"
df_phy_run4 <- merge.all[,c("type","phy_run4")]
df_phy_run4$class <- "run"
colnames(df_phy_run4)[2] <- "phy"

data_frames_to_merge <- list(
  df_phy_loop1,
  df_phy_loop2,
  df_phy_loop3,
  df_phy_run1,
  df_phy_run2,
  df_phy_run3,
  df_phy_run4
)

# 使用 Reduce 函数按列合并数据框
df <- Reduce(rbind, data_frames_to_merge)
df$type <- factor(df$type,levels =c("non_eG4","kc_specific","s2_specific","overlap","merge"))

# 查看表格中是否有缺失值
if (anyNA(df)) {
  cat("表格中存在缺失值。\n")
} else {
  cat("表格中没有缺失值。\n")
}

ggplot(data = df,aes(x=type,y=phy,fill=class)) +
  geom_boxplot(notch = TRUE,outlier.colour = "white") +
  scale_fill_manual(values = c("#74ADD1","#D6604D")) +
  cowplot::theme_half_open()+
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank()) +
  ylab("PhyloP score") +
  stat_compare_means(aes(group = class),method = "t.test",label = "p.signif",label.y = 2.5) +
  coord_cartesian(ylim = c(-1.6,2.6)) 
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"G4LoopRun.PhyloP.pdf"),
       device = "pdf",width = 7,height = 5)
