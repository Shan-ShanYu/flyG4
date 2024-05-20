#cat /home/qians/G4/Drosophila/G4FS/F-run/H4K16ac_MSL/male_W_peaks.bed /home/qians/G4/Drosophila/G4FS/F-run/H4K16ac_MSL/female_W_peaks.bed | bedtools sort -i - |bedtools merge -i - > 006.1.H4K16ac.bothsex.merge.bed

# bedtools intersect -a 006.1.genome.1000.window -b "/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.kc_specific.chr.bed" -wa -c | bedtools intersect -a - -b "/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.s2_specific.chr.bed" -wa -c | bedtools intersect -a - -b "/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.overlap.chr.bed" -wa -c | bedtools intersect -a - -b "/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.merge.chr.bed" -wa -c > 006.1.genome1000.eG4.bed
# bedtools intersect -a 006.1.HAS.1000.bed -b "/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.kc_specific.chr.bed" -wa -c | bedtools intersect -a - -b "/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.s2_specific.chr.bed" -wa -c | bedtools intersect -a - -b "/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.overlap.chr.bed" -wa -c | bedtools intersect -a - -b "/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.merge.chr.bed" -wa -c > 006.1.HAS1000.eG4.bed
# bedtools intersect -a 006.1.H4K16ac.bothsex.merge.1000.bed -b "/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.kc_specific.chr.bed" -wa -c | bedtools intersect -a - -b "/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.s2_specific.chr.bed" -wa -c | bedtools intersect -a - -b "/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.overlap.chr.bed" -wa -c | bedtools intersect -a - -b "/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.merge.chr.bed" -wa -c > 006.1.H4K16ac.bothsex1000.eG4.bed

#### H4k16ac 在所有染色体上 ####
rm(list = ls());gc();rm(list = ls())
Num = "006.1."
window = fread("/home/yuss/flyG4/result/H4k16ac/006.1.genome1000.eG4.bed") %>% as.data.frame()
table(window$V1)
# chr2L chr2R chr3L chr3R  chr4  chrX  chrY 
# 23514 25287 28111 32080  1349 23543  3668 
subset_window <- subset(window, window$V1 == "chr2L")
total_sum <- sum(subset_window$V4)
558/23514 #G4数量除以window数，也就是得到了平均每个window，G4的数量
tapply(window$V4,window$V1,mean) #kc
# chr2L       chr2R       chr3L       chr3R        chr4        chrX        chrY 
# 0.023730544 0.025546724 0.023264914 0.032169576 0.000000000 0.068470458 0.001090513 
tapply(window$V5,window$V1,mean) #s2
# chr2L       chr2R       chr3L       chr3R        chr4        chrX        chrY 
# 0.027217828 0.018744810 0.022589022 0.015305486 0.005930319 0.010194113 0.000000000 
tapply(window$V6,window$V1,mean) #overlap
# chr2L       chr2R       chr3L       chr3R        chr4        chrX        chrY 
# 0.070298546 0.090164907 0.073956814 0.089588529 0.005189029 0.042900225 0.000000000 
tapply(window$V7,window$V1,mean) #merge
# chr2L       chr2R       chr3L       chr3R        chr4        chrX        chrY 
# 0.121246917 0.134456440 0.119810750 0.137063591 0.011119348 0.121564796 0.001090513 
mean(window[window$V1!="chr4"&window$V1!="chrX"&window$V1!="chrY","V7"]) #eG4在常染色体上平均密度是1.251628
# 0.1285966
tapply(window$V7,window$V1,sum)
table(window$V1)
# mean(window$V7[c(which(grepl("R",window$V1)),which(grepl("L",window$V1)))])
#* KC167--------------------------------------------------------------------------------------------------
has = fread("/home/yuss/flyG4/result/H4k16ac/006.1.HAS1000.eG4.bed") %>% as.data.frame()
sum(has$V4)/nrow(has) #kc 0.01167315 #G4数量除以window数
nrow(has) #257
p = matrix(0,nrow = 5,ncol = 3)
for (condition in 1) {
  for (i in 1:5) {
    chr = paste0("chr",c("2L","2R","3L","3R","X"))[i]
    p[i,condition] = wilcox.test(window[window$V1==chr,condition+3],has[,condition+3])$p.value
  }
}
#kc 0.002028843

h4k16ac <- fread("/home/yuss/flyG4/result/H4k16ac/006.1.H4K16ac.bothsex1000.eG4.bed") %>% as.data.frame()
h4k16ac$len <- h4k16ac$V3-h4k16ac$V2
sum(h4k16ac$V4)/sum(h4k16ac$len)*10^3 #kc 0.04030094 计算G4的密度*10^3，表示1000bp上G4的数量，即平均每个window G4的数量
p = matrix(0,nrow = 5,ncol = 3)
for (condition in 1) {
  for (i in 1:5) {
    chr = paste0("chr",c("2L","2R","3L","3R","X"))[i]
    p[i,condition] =wilcox.test(window[window$V1==chr,condition+3],h4k16ac[,condition+3])$p.value
  }
}
#kc 1.700170e-57
df_kc = data.frame(name=c("2L","2R","3L","3R","X","HAS","H4K16ac"),
                num=c(as.numeric(tapply(window$V4,window$V1,mean))[c(1:4,6)],sum(has$V4)/nrow(has),sum(h4k16ac$V4)/sum(h4k16ac$len)*10^3)) 
#* S2--------------------------------------------------------------------------------------------------
sum(has$V5)/nrow(has) #s2 0.0233463
nrow(has) #257
p = matrix(0,nrow = 5,ncol = 3)
for (condition in 1) {
  for (i in 1:5) {
    chr = paste0("chr",c("2L","2R","3L","3R","X"))[i]
    p[i,condition] = wilcox.test(window[window$V1==chr,condition+4],has[,condition+4])$p.value
  }
}
#s2 0.2342594

sum(h4k16ac$V5)/sum(h4k16ac$len)*10^3 #s2 0.01621214
p = matrix(0,nrow = 5,ncol = 3)
for (condition in 1) {
  for (i in 1:5) {
    chr = paste0("chr",c("2L","2R","3L","3R","X"))[i]
    p[i,condition] =wilcox.test(window[window$V1==chr,condition+4],h4k16ac[,condition+4])$p.value
  }
}
#s2 1.037506e-04
df_s2 = data.frame(name=c("2L","2R","3L","3R","X","HAS","H4K16ac"),
                   num=c(as.numeric(tapply(window$V5,window$V1,mean))[c(1:4,6)],sum(has$V5)/nrow(has),sum(h4k16ac$V5)/sum(h4k16ac$len)*10^3))
#* overlap-----------------------------------------------------------------------------------------------
sum(has$V6)/nrow(has) #overlap 0.03501946
nrow(has) #257
p = matrix(0,nrow = 5,ncol = 3)
for (condition in 1) {
  for (i in 1:5) {
    chr = paste0("chr",c("2L","2R","3L","3R","X"))[i]
    p[i,condition] = wilcox.test(window[window$V1==chr,condition+5],has[,condition+5])$p.value
  }
}
#overlap 0.698717662

sum(h4k16ac$V6)/sum(h4k16ac$len)*10^3 #overlap 0.06924574
p = matrix(0,nrow = 5,ncol = 3)
for (condition in 1) {
  for (i in 1:5) {
    chr = paste0("chr",c("2L","2R","3L","3R","X"))[i]
    p[i,condition] =wilcox.test(window[window$V1==chr,condition+5],h4k16ac[,condition+5])$p.value
  }
} ##实际上是在比较 window 数据框中与特定染色体（如 chr2L、chr2R 等）相关的子集与 h4k16ac 数据框的第六列之间是否存在显著差异。
#overlap 1.584152e-16
df_overlap = data.frame(name=c("2L","2R","3L","3R","X","HAS","H4K16ac"),
                   num=c(as.numeric(tapply(window$V6,window$V1,mean))[c(1:4,6)],sum(has$V6)/nrow(has),sum(h4k16ac$V6)/sum(h4k16ac$len)*10^3))

#* merge-----------------------------------------------------------------------------------------------
sum(has$V7)/nrow(has) #merge 0.07003891
p = matrix(0,nrow = 5,ncol = 3)
for (condition in 1) {
  for (i in 1:5) {
    chr = paste0("chr",c("2L","2R","3L","3R","X"))[i]
    p[i,condition] = wilcox.test(window[window$V1==chr,condition+6],has[,condition+6])$p.value
  }
}
#merge 0.039758793

sum(h4k16ac$V7)/sum(h4k16ac$len)*10^3 #merge 0.1257588
p = matrix(0,nrow = 5,ncol = 3)
for (condition in 1) {
  for (i in 1:5) {
    chr = paste0("chr",c("2L","2R","3L","3R","X"))[i]
    p[i,condition] =wilcox.test(window[window$V1==chr,condition+6],h4k16ac[,condition+6])$p.value
  }
} ##实际上是在比较 window 数据框中与特定染色体（如 chr2L、chr2R 等）相关的子集与 h4k16ac 数据框的第六列之间是否存在显著差异。
#merge 1.022399e-03
df_merge = data.frame(name=c("2L","2R","3L","3R","X","HAS","H4K16ac"),
                        num=c(as.numeric(tapply(window$V7,window$V1,mean))[c(1:4,6)],sum(has$V7)/nrow(has),sum(h4k16ac$V7)/sum(h4k16ac$len)*10^3))



#### H4k16ac 在X染色体上 ####
rm(list = ls());gc();rm(list = ls())
Num = "006.1."
window = fread("/home/yuss/flyG4/result/H4k16ac/006.1.genome1000.eG4.bed") %>% as.data.frame()
table(window$V1)
# chr2L chr2R chr3L chr3R  chr4  chrX  chrY 
# 23514 25287 28111 32080  1349 23543  3668 
#* KC167--------------------------------------------------------------------------------------------------
has = fread("/home/yuss/flyG4/result/H4k16ac/006.1.HAS1000.eG4.bed") %>% as.data.frame()
sum(has$V4)/nrow(has) #kc 0.01167315
nrow(has) #257
p = matrix(0,nrow = 5,ncol = 3)
for (condition in 1) {
  for (i in 1:5) {
    chr = paste0("chr",c("2L","2R","3L","3R","X"))[i]
    p[i,condition] = wilcox.test(window[window$V1==chr,condition+3],has[,condition+3])$p.value
  }
}
#kc 0.002028843

h4k16ac <- fread("/home/yuss/flyG4/result/H4k16ac/006.1.H4K16ac.bothsex1000.eG4.bed") %>% as.data.frame()
h4k_chrx <- h4k16ac[h4k16ac$V1=="chrX",]
h4k_chrx$len <- h4k_chrx$V3-h4k_chrx$V2
sum(h4k_chrx$V4)/sum(h4k_chrx$len)*10^3 #kc 0.06133505 计算G4的密度*10^3，表示1000bp上G4的数量
p = matrix(0,nrow = 5,ncol = 3)
for (condition in 1) {
  for (i in 1:5) {
    chr = paste0("chr",c("2L","2R","3L","3R","X"))[i]
    p[i,condition] =wilcox.test(window[window$V1==chr,condition+3],h4k_chrx[,condition+3])$p.value
  }
}
#kc 3.232585e-04
df_kc_x = data.frame(name=c("2L","2R","3L","3R","X","HAS","h4k_chrx"),
                   num=c(as.numeric(tapply(window$V4,window$V1,mean))[c(1:4,6)],sum(has$V4)/nrow(has),sum(h4k_chrx$V4)/sum(h4k_chrx$len)*10^3)) 

#* S2--------------------------------------------------------------------------------------------------
sum(has$V5)/nrow(has) #s2 0.0233463
nrow(has) #257
p = matrix(0,nrow = 5,ncol = 3)
for (condition in 1) {
  for (i in 1:5) {
    chr = paste0("chr",c("2L","2R","3L","3R","X"))[i]
    p[i,condition] = wilcox.test(window[window$V1==chr,condition+4],has[,condition+4])$p.value
  }
}
#s2 0.2342594

sum(h4k_chrx$V5)/sum(h4k_chrx$len)*10^3 #s2 0.01169029
p = matrix(0,nrow = 5,ncol = 3)
for (condition in 1) {
  for (i in 1:5) {
    chr = paste0("chr",c("2L","2R","3L","3R","X"))[i]
    p[i,condition] =wilcox.test(window[window$V1==chr,condition+4],h4k_chrx[,condition+4])$p.value
  }
}
#s2 5.729266e-01
df_s2_x = data.frame(name=c("2L","2R","3L","3R","X","HAS","h4k_chrx"),
                   num=c(as.numeric(tapply(window$V5,window$V1,mean))[c(1:4,6)],sum(has$V5)/nrow(has),sum(h4k_chrx$V5)/sum(h4k_chrx$len)*10^3))
#* overlap-----------------------------------------------------------------------------------------------
sum(has$V6)/nrow(has) #overlap 0.03501946
nrow(has) #257
p = matrix(0,nrow = 5,ncol = 3)
for (condition in 1) {
  for (i in 1:5) {
    chr = paste0("chr",c("2L","2R","3L","3R","X"))[i]
    p[i,condition] = wilcox.test(window[window$V1==chr,condition+5],has[,condition+5])$p.value
  }
}
#overlap 0.698717662

sum(h4k_chrx$V6)/sum(h4k_chrx$len)*10^3 #overlap 0.04520245
p = matrix(0,nrow = 5,ncol = 3)
for (condition in 1) {
  for (i in 1:5) {
    chr = paste0("chr",c("2L","2R","3L","3R","X"))[i]
    p[i,condition] =wilcox.test(window[window$V1==chr,condition+5],h4k_chrx[,condition+5])$p.value
  }
} ##实际上是在比较 window 数据框中与特定染色体（如 chr2L、chr2R 等）相关的子集与 h4k_chrx 数据框的第六列之间是否存在显著差异。
#overlap 8.740450e-01
df_overlap_x = data.frame(name=c("2L","2R","3L","3R","X","HAS","h4k_chrx"),
                        num=c(as.numeric(tapply(window$V6,window$V1,mean))[c(1:4,6)],sum(has$V6)/nrow(has),sum(h4k_chrx$V6)/sum(h4k_chrx$len)*10^3))

#* merge-----------------------------------------------------------------------------------------------
sum(has$V7)/nrow(has) #merge 0.07003891
p = matrix(0,nrow = 5,ncol = 3)
for (condition in 1) {
  for (i in 1:5) {
    chr = paste0("chr",c("2L","2R","3L","3R","X"))[i]
    p[i,condition] = wilcox.test(window[window$V1==chr,condition+6],has[,condition+6])$p.value
  }
}
#merge 0.039758793

sum(h4k_chrx$V7)/sum(h4k_chrx$len)*10^3 #merge 0.1182278
p = matrix(0,nrow = 5,ncol = 3)
for (condition in 1) {
  for (i in 1:5) {
    chr = paste0("chr",c("2L","2R","3L","3R","X"))[i]
    p[i,condition] =wilcox.test(window[window$V1==chr,condition+6],h4k_chrx[,condition+6])$p.value
  }
} ##实际上是在比较 window 数据框中与特定染色体（如 chr2L、chr2R 等）相关的子集与 h4k_chrx 数据框的第六列之间是否存在显著差异。
#merge 2.449976e-02
df_merge_x = data.frame(name=c("2L","2R","3L","3R","X","HAS","h4k_chrx"),
                      num=c(as.numeric(tapply(window$V7,window$V1,mean))[c(1:4,6)],sum(has$V7)/nrow(has),sum(h4k_chrx$V7)/sum(h4k_chrx$len)*10^3))


# 比较了师兄G4seq的文件和自己的文件，两者是一样的，不过自己的文件是排序后的
# sort /home/qians/Quadruplex/Input/Data/G4seq/dm6.K.bed | less
# less "/home/yuss/flyG4/data/G4-seq/GSM3003541_Drosophila_all_w15_th-1.K.bed"
# 
# #H4K16ac
# ##1.复制female male的peaks.bed文件
# (base) yuss@ubuntu:~/flyG4/result/H4k16ac$ cp /home/qians/G4/Drosophila/G4FS/F-run/H4K16ac_MSL/female_W_peaks.bed 006.2.female_W_peaks.bed
# (base) yuss@ubuntu:~/flyG4/result/H4k16ac$ cp /home/qians/G4/Drosophila/G4FS/F-run/H4K16ac_MSL/male_W_peaks.bed 006.2.male_W_peaks.bed
# 
# ##2.female male的peaks.bed文件与G4取交集，得到female male的H4K16ac文件
# for i in male female; do bedtools intersect -a /home/yuss/flyG4/result/H4k16ac/006.2.${i}_W_peaks.bed -b "/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.merge.chr.bed" -wa -c | bedtools  intersect -a - -b "/home/yuss/flyG4/data/G4-seq/GSM3003541_Drosophila_all_w15_th-1.K.bed" -wa -c | bedtools intersect -a - -b "/home/yuss/flyG4/data/G4-seq/GSM3003542_Drosophila_all_w15_th-1.PDS.bed" -wa -c > 006.2.${i}.H4K16ac.MergedG4.K.PDS.bed; done
# 
# ##3.female male的Flank1kb文件
# for i in male female; do bedtools flank -i /home/yuss/flyG4/result/H4k16ac/006.2.${i}_W_peaks.bed -g "/home/yuss/flyG4/data/ref/dm6.genome" -b 1000 | bedtools intersect -a - -b "/home/yuss/flyG4/result/PQS/5.type.bed/chr.bed/001.2.merge.chr.bed" -wa -c | bedtools  intersect -a - -b "/home/yuss/flyG4/data/G4-seq/GSM3003541_Drosophila_all_w15_th-1.K.bed" -wa -c | bedtools intersect -a - -b "/home/yuss/flyG4/data/G4-seq/GSM3003542_Drosophila_all_w15_th-1.PDS.bed" -wa -c > 006.2.${i}.Flank1kb.MergedG4.K.PDS.bed; done

rm(list = ls());gc();rm(list = ls())
Num = "006.2."
#### H4k16ac ####
#*MergeG4----------------------------------------------------------------------------------
df = c()
for (sex in c("male","female")) {
  for (region in c("H4K16ac","Flank1kb")) {
    a = fread(paste0("/home/yuss/flyG4/result/H4k16ac/",Num,sex,".",region,".MergedG4.K.PDS.bed")) %>% as.data.frame()
    a %<>% dplyr::filter(.,V1 %in% paste0("chr",c("2L","2R","3L","3R","X")))
    a$V1 %<>% gsub("chr","",.)
    a[a$V1!="X","V1"]="A"
    a$type = paste(region,sex,a$V1,sep = "_")
    df = rbind(df,a)
  }
}
df$len = df$V3-df$V2

den = tapply(df$V10, df$type, sum)/tapply(df$len, df$type, sum) %>% as.data.frame()
den$. =den$. * 10^3
den$type = row.names(den) %>% gsub("_female_","_F_",.) %>% gsub("_male_","_M_",.) %>% gsub("1kb","",.)
den$type %<>% factor(.,levels = c(paste(rep("H4K16ac",4),rep(c("F","M"),each=2),rep(c("A","X"),2),sep="_"),
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

#### HAS ####
# 在HAS位点 G4 density （HAS是在X染色体上的，并且剂量补偿是在雄性果蝇中发生的，只有一个性别和一种染色体，无法与常染色体或者雌性相比，做不了）

#### H4k16ac 在X染色体上 雄性####
window = fread("/home/yuss/flyG4/result/H4k16ac/006.1.genome1000.eG4.bed") %>% as.data.frame()
h4k16ac.male <- fread("/home/yuss/flyG4/result/H4k16ac/006.1.H4K16ac.male1000.eG4.bed") %>% as.data.frame()
h4k16ac.male$len <- h4k16ac.male$V3-h4k16ac.male$V2
h4kmale_chrx <- h4k16ac.male[h4k16ac.male$V1=="chrX",]
nrow(h4kmale_chrx)
sum(h4kmale_chrx$V4)/sum(h4kmale_chrx$len)*10^3 ##0.1165653 计算G4的密度*10^3，表示1000bp上G4的数量
p = matrix(0,nrow = 5,ncol = 3)
for (condition in 1) {
  for (i in 1:5) {
    chr = paste0("chr",c("2L","2R","3L","3R","X"))[i]
    p[i,condition] =wilcox.test(window[window$V1==chr,condition+3],h4kmale_chrx[,condition+3])$p.value
  }
}
#male 8.110344e-34

#### H4k16ac 在X染色体上 雌性####
h4k16ac.female <- fread("/home/yuss/flyG4/result/H4k16ac/006.1.H4K16ac.female1000.eG4.bed") %>% as.data.frame()
h4k16ac.female$len <- h4k16ac.female$V3-h4k16ac.female$V2
h4kfemale_chrx <- h4k16ac.female[h4k16ac.female$V1=="chrX",]
nrow(h4kfemale_chrx)
sum(h4kfemale_chrx$V4)/sum(h4kfemale_chrx$len)*10^3 ##0.1387177
p = matrix(0,nrow = 5,ncol = 3)
for (condition in 1) {
  for (i in 1:5) {
    chr = paste0("chr",c("2L","2R","3L","3R","X"))[i]
    p[i,condition] =wilcox.test(window[window$V1==chr,condition+3],h4kfemale_chrx[,condition+3])$p.value
  }
}
#female 1.042533e-17

#### kc_all/s2_all G4 ####
rm(list = ls());gc();rm(list = ls())
Num = "006.1."
window = fread("/home/yuss/flyG4/result/H4k16ac/006.1.genome1000.kc_all.s2_all.eG4.bed") %>% as.data.frame()
table(window$V1)
# chr2L chr2R chr3L chr3R  chr4  chrX  chrY 
# 23514 25287 28111 32080  1349 23543  3668 
tapply(window$V4,window$V1,mean) #kcall
# chr2L       chr2R       chr3L       chr3R        chr4        chrX        chrY 
# 0.094029089 0.115711630 0.097221728 0.121758105 0.005189029 0.111370683 0.001090513
tapply(window$V5,window$V1,mean) #s2all
# chr2L      chr2R      chr3L      chr3R       chr4       chrX       chrY 
# 0.09751637 0.10890972 0.09654584 0.10489401 0.01111935 0.05309434 0.00000000 

window$chr1 <- ifelse(window$V1=="chrX","X","A")
table(window$chr1)
# A      X 
# 114009  23543 
tapply(window$V4,window$chr1,mean) #kcall
# A         X 
# 0.1033866 0.1113707
tapply(window$V5,window$chr1,mean) #s2all
# A          X 
# 0.09772036 0.05309434 
sum(window[window$V1 == "chrX", "V5"])
window$len <- window$V3-window$V2
sum(window[window$V1 == "chrX", "len"])

#### A X ####
A <- window[window$chr1=="A",]
X <- window[window$chr1=="X",]
# 初始化 p 矩阵
p <- matrix(0, nrow = 5, ncol = 2)
# 进行 Wilcoxon 秩和检验
p <- wilcox.test(A[,4], X[,4])$p.value ###kcall 0.0453
p <- wilcox.test(A[,5], X[,5])$p.value ###s2all 8.927926e-87

#### has s2 ####
has = fread("/home/yuss/flyG4/result/H4k16ac/006.1.HAS1000.s2_all.eG4.bed") %>% as.data.frame()
sum(has$V4)/nrow(has) #s2 0.05836576
nrow(has) #257
p = matrix(0,nrow = 5,ncol = 3)
for (condition in 1) {
  for (i in 1:5) {
    chr = paste0("chr",c("2L","2R","3L","3R","X"))[i]
    p[i,condition] = wilcox.test(window[window$V1==chr,condition+3],has[,condition+3])$p.value
  }
}
#s2 0.01899724

#### H4k16ac雄性与s2交集 在X染色体上 ####
h4k16ac.male <- fread("/home/yuss/flyG4/result/H4k16ac/006.1.H4K16ac.male1000.s2_all.bed") %>% as.data.frame()
h4k16ac.male$len <- h4k16ac.male$V3-h4k16ac.male$V2
h4kmale_chrx <- h4k16ac.male[h4k16ac.male$V1=="chrX",]
sum(h4kmale_chrx$V4)/sum(h4kmale_chrx$len)*10^3 ##0.05466286 计算G4的密度*10^3，表示1000bp上G4的数量
p = matrix(0,nrow = 5,ncol = 3)
for (condition in 1) {
  for (i in 1:5) {
    chr = paste0("chr",c("2L","2R","3L","3R","X"))[i]
    p[i,condition] =wilcox.test(window[window$V1==chr,condition+3],h4kmale_chrx[,condition+3])$p.value
  }
}
#male 3.619614e-61

#### H4k16ac雌性与kc交集 在X染色体上 ####
h4k16ac.female <- fread("/home/yuss/flyG4/result/H4k16ac/006.1.H4K16ac.female1000.kc_all.bed") %>% as.data.frame()
h4k16ac.female$len <- h4k16ac.female$V3-h4k16ac.female$V2
h4kfemale_chrx <- h4k16ac.female[h4k16ac.female$V1=="chrX",]
sum(h4kfemale_chrx$V4)/sum(h4kfemale_chrx$len)*10^3 ##0.126107
p = matrix(0,nrow = 5,ncol = 3)
for (condition in 1) {
  for (i in 1:5) {
    chr = paste0("chr",c("2L","2R","3L","3R","X"))[i]
    p[i,condition] =wilcox.test(window[window$V1==chr,condition+3],h4kfemale_chrx[,condition+3])$p.value
  }
}
#female 0.1840751815


#### 画图 ####
##kc画图
mean(window[window$V1!="chr4"&window$V1!="chrX"&window$V1!="chrY","V4"]) #eG4在常染色体上平均密度是0.1080446
##chx 0.111370683
##h4kfemale X 0.126107
##h4kfemale A 0.1163629
h4kfemale_chrA <- h4k16ac.female[h4k16ac.female$V1!="chrX"&h4k16ac.female$V1!="chr4"&h4k16ac.female$V1!="chrY",]
sum(h4kfemale_chrA$V4)/sum(h4kfemale_chrA$len)*10^3 ##0.1163629


##s2画图
mean(window[window$V1!="chr4"&window$V1!="chrX"&window$V1!="chrY","V5"]) #eG4在常染色体上平均密度是0.1020809
##chx 0.05309434
##h4kmale X 0.05466286
##h4kmale A 0.1183512
h4kmale_chrA <- h4k16ac.male[h4k16ac.male$V1!="chrX"&h4k16ac.male$V1!="chr4"&h4k16ac.male$V1!="chrY",]
sum(h4kmale_chrA$V4)/sum(h4kmale_chrA$len)*10^3 ##0.1183512 计算G4的密度*10^3，表示1000bp上G4的数量
