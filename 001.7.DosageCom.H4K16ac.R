# #H4K16ac
# ##复制female male的文件
# for i in male female; do bedtools intersect -a /home/yuss/flyG4/result/PQS/001.7.${i}_W_peaks.bed -b "/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.merge.chr.bed"  -wa -c > 001.7.${i}.H4K16ac.MergedG4.bed; done
# #Flank 1 kb
# ##复制基因组大小文件
# cp /home/qians/Quadruplex/Input/Ref/Fly/dm6.genome /home/yuss/flyG4/data/ref/
#   vi /home/qians/Quadruplex/Input/Ref/Fly/dm6.genome ##加上chr
# for i in male female; do bedtools flank -i /home/yuss/flyG4/result/PQS/001.7.${i}_W_peaks.bed -g "/home/yuss/flyG4/data/ref/dm6.genome" -b 1000 | bedtools intersect -a - -b "/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.merge.chr.bed" -wa -c > 001.7.${i}.Flank1kb.MergedG4.bed; done
rm(list = ls());gc();rm(list = ls())
Num = "001.7."

#*MergedG4---------------------------------------------------------
Num = "001.7."
df = c()
for (sex in c("male","female")) {
  for (region in c("H4K16ac","Flank1kb")) {
    a = fread(paste0("/home/yuss/flyG4/result/PQS/",Num,sex,".",region,".MergedG4.K.PDS.bed")) %>% as.data.frame()
    a %<>% dplyr::filter(.,V1 %in% paste0("chr",c("2L","2R","3L","3R","X")))
    a$V1 %<>% gsub("chr","",.)
    a[a$V1!="X","V1"]="A"
    a$type = paste(region,sex,a$V1,sep = "_")
    df = rbind(df,a)
  }
}
df$len = df$V3 - df$V2

den = tapply(df$V10, df$type, sum)/tapply(df$len, df$type, sum) %>% as.data.frame()
den$. = den$. * 10^3
den$type = row.names(den) %>% gsub("_female_","_F_",.) %>% gsub("_male_","_M_",.) %>% gsub("1kb","",.)
den$type %<>% factor(.,levels = c(paste(rep("H4K16ac",4),rep(c(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAWElEQVR42mNgGPTAxsZmJsVqQApgmGw1yApwKcQiT7phRBuCzzCSDSHGMKINIeDNmWQlA2IigKJwIssQkHdINgxfmBBtGDEBS3KCxBc7pMQgMYE5c/AXPwAwSX4lV3pTWwAAAABJRU5ErkJggg=="F","M"),each=2),rep(c("A","X"),2),sep="_"),
                                  paste(rep("Flank",4),rep(c("F","M"),each=2),rep(c("A","X"),2),sep="_")))
ggplot(den,aes(type,.))+geom_bar(stat = "identity",position = "dodge")+
  theme_bw()+ ylab("Density (No./kb)")+#+xlab("Group")+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=10,colour = "black",angle = 45,hjust = 1,vjust = 1),axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,colour = "black"),axis.title.y = element_text(size=12))+
  scale_y_continuous(expand = c(0,0),limits = c(0,0.157))+
  geom_hline(yintercept = mean(den$.),color="#d6604d",size=1.2,linetype="dotted")+
  geom_vline(xintercept = 4.5,color="#595959",size=0.8) +
  geom_signif(y_position=0.143, xmin=3, xmax=4,
              annotation=c("***"), tip_length=0.06)#p-value < 2.2e-16
  # geom_text(data=den[den$sig,],aes(label = "*"), vjust = -0.5)
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"DosageCom.DensityH4Flank.MergedG4.pdf"),
       device = "pdf",width = 4.2,height = 3)

##不画flank区域
tapply(window$V7,window$V1,mean) #merge
# chr2L       chr2R       chr3L       chr3R        chr4        chrX        chrY 
# 0.121246917 0.134456440 0.119810750 0.137063591 0.011119348 0.121564796 0.001090513 
mean(window[window$V1!="chr4"&window$V1!="chrX"&window$V1!="chrY","V7"]) #eG4在常染色体上平均密度是1.251628
# 0.1285966

den1 <- den[5:8,]
#常染色体、X染色体和HAS的密度
new_values <- data.frame(
  "." = c(0.1285966, 0.121564796, 0.07003891),
  type = c("Autosome", "X chromosome", "HAS")
)

# 将新的数值添加到原有表格下面
den1_extended <- rbind(den1, new_values)
den1_extended$type <- factor(den1_extended$type,levels=c("Autosome","X chromosome","HAS","H4K16ac_F_A","H4K16ac_F_X","H4K16ac_M_A","H4K16ac_M_X"))
# 打印新的表格
print(den1_extended)

ggplot(den1_extended,aes(type,.))+geom_bar(stat = "identity",position = "dodge")+
  theme_bw()+ ylab("Density (No./kb)")+#+xlab("Group")+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=10,colour = "black",angle = 45,hjust = 1,vjust = 1),axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,colour = "black"),axis.title.y = element_text(size=12))+
  scale_y_continuous(expand = c(0,0),limits = c(0,0.17))+
  # geom_hline(yintercept = mean(den$.),color="#d6604d",size=1.2,linetype="dotted")+
  geom_vline(xintercept = 3.5,color="#595959",size=0.5,linetype="dotted") +
  geom_signif(y_position=0.142, xmin=2, xmax=5,
              annotation=c("****"), tip_length=0.06) + #p-value < 2.2e-16
  geom_signif(y_position=0.155, xmin=2, xmax=7,
              annotation=c("****"), tip_length=0.06) 

# geom_text(data=den[den$sig,],aes(label = "*"), vjust = -0.5)
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"DosageCom.DensityH4.MergedG4.pdf"),
       device = "pdf",width = 5,height = 3.5)


# ##密度
# # H4K16ac_male_A
# # 1.350894e-04
# # H4K16ac_male_X
# # 1.134742e-04
# #sum density 24856365
# tapply(df$len, df$type, sum)
# # H4K16ac_male_A 5914601
# # H4K16ac_male_X 12293540
# #18208141
# data <- matrix(c(13508942,24856365*5914601/18208141,11347423,24856365*12293540/18208141),nrow = 2)
# colnames(data) <- c("A","X")
# rownames(data) <- c("den","predict_den")
# fisher.test(data) # p-value < 2.2e-16
# 
# ##数量
# tapply(df$V10, df$type, sum)
# # H4K16ac_male_A 799
# # H4K16ac_male_X 1395
# #sum num 2194
# data1 <- matrix(round(c(799,2194*5914601/18208141,1395,2194*12293540/18208141),0),nrow = 2)               
# colnames(data1) <- c("A","X")
# rownames(data1) <- c("num","predict_num")
# fisher.test(data1) #0.006922

##修改后正确的是:假设在X和常染色体，G4的密度一样
# H4K16ac_male_A 799
# H4K16ac_male_X 1395
#sum num 2194
#sum length 18208141
data <- matrix(round((c(1.350894e-04,2194/18208141,1.134742e-04,2194/18208141)*10000000),0),nrow = 2)
colnames(data) <- c("A","X")
rownames(data) <- c("den","predict_den")
fisher.test(data) # p-value < 2.2e-16
#*K-----------------------------------------------------------------------------
den = tapply(df$V11, df$type, sum)/tapply(df$len, df$type, sum) %>% as.data.frame()
den$. = den$. * 10^3
den$type = row.names(den) %>% gsub("_female_","_F_",.) %>% gsub("_male_","_M_",.) %>% gsub("1kb","",.)
den$type %<>% factor(.,levels = c(paste(rep("H4K16ac",4),rep(c("F","M"),each=2),rep(c("A","X"),2),sep="_"),
                                  paste(rep("Flank",4),rep(c("F","M"),each=2),rep(c("A","X"),2),sep="_")))
ggplot(den,aes(type,.))+geom_bar(stat = "identity",position = "dodge")+
  theme_bw()+ ylab("Density (No./kb)")+#+xlab("Group")+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=10,colour = "black",angle = 45,hjust = 1,vjust = 1),axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,colour = "black"),axis.title.y = element_text(size=12))+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(den$.)*1.1))+
  geom_hline(yintercept = mean(den$.),color="#d6604d",size=1.2,linetype="dotted")+
  geom_vline(xintercept = 4.5,color="#595959",size=0.8)
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"DosageCom.DensityH4Flank.G4seqK.pdf"),
       device = "pdf",width = 4.2,height = 3)

#*PDS-----------------------------------------
den = tapply(df$V12, df$type, sum)/tapply(df$len, df$type, sum) %>% as.data.frame()
den$. = den$. * 10^3
den$type = row.names(den) %>% gsub("_female_","_F_",.) %>% gsub("_male_","_M_",.) %>% gsub("1kb","",.)
den$type %<>% factor(.,levels = c(paste(rep("H4K16ac",4),rep(c("F","M"),each=2),rep(c("A","X"),2),sep="_"),
                                  paste(rep("Flank",4),rep(c("F","M"),each=2),rep(c("A","X"),2),sep="_")))
ggplot(den,aes(type,.))+geom_bar(stat = "identity",position = "dodge")+
  theme_bw()+ ylab("Density (No./kb)")+#+xlab("Group")+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=10,colour = "black",angle = 45,hjust = 1,vjust = 1),axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,colour = "black"),axis.title.y = element_text(size=12))+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(den$.)*1.1))+
  geom_hline(yintercept = mean(den$.),color="#d6604d",size=1.2,linetype="dotted")+
  geom_vline(xintercept = 4.5,color="#595959",size=0.8)
ggsave(filename = paste0("/home/yuss/flyG4/result/PQS/Picture/",Num,"DosageCom.DensityH4Flank.G4seqPDS.pdf"),
       device = "pdf",width = 4.2,height = 3)

#*MergedG4-----------------------------------------
