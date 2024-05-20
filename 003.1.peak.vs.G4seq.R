rm(list = ls());gc();rm(list = ls())#清空
Num = "003.1."
#### PQS与peak,G4seq取交集后表的整合 ####
pqs.peak <- fread("/home/yuss/flyG4/result/data.reliability/003.1.pqs.peak.bed")
pqs.G4seq <- fread("/home/yuss/flyG4/result/data.reliability/003.1.pqs.G4seq.bed")
colnames(pqs.peak)[7] <- 'peak'
colnames(pqs.G4seq)[7] <- 'G4seq'
df <- bind_cols(pqs.peak,pqs.G4seq$G4seq)
colnames(df) <- c("chr","start","end","pqs.id","score","strand","peak","G4seq")
df$peak <- ifelse(df$peak==0,0,1)
df$G4seq <- if_else(df$G4seq==0,0,2)
df$sum <- rowSums(df[,7:8])
peak.specific <- df[df$sum==1,]
G4seq.specific <- df[df$sum==2,]
overlap <- df[df$sum==3,]
non_eG4 <- df[df$sum==0,]
peak.specific$type <- "peak.specific"
G4seq.specific$type <- "G4seq.specific"
overlap$type <- "overlap"
non_eG4$type <- "non_eG4"
df.intersect <- bind_rows(non_eG4,peak.specific,G4seq.specific,overlap)
write.table(df.intersect,file = paste0("/home/yuss/flyG4/result/data.reliability/",Num,"df.intersect.bed"),
            sep = '\t', col.names = T,row.names = F,quote = F)
table(df.intersect$type)
##画成一个表格
data <- data.frame(
  PQS_intersect_with_peak = c(17072,3907),
  PQS_not_intersect_with_peak =c(17880,5058))
rownames(data) <- c("PQS_intersect_with_G4seq","PQS_not_intersect_with_G4seq")

#### 结构稳定性 ####
df.intersect$type <- factor(df.intersect$type,levels = c("non_eG4","peak.specific","G4seq.specific","overlap"))
median.score <- aggregate(score ~ type, df.intersect, median)
library(ggpubr)
ggplot(data = df.intersect,aes(x=type,y=score,fill=type))+
  geom_boxplot(notch = TRUE,outlier.colour = "white")
