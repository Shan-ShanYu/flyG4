rm(list = ls());gc();rm(list = ls())

#*KC S2合并的G4做-----------------------------------------------------------------------------------------
window = fread("/home/yuss/flyG4/result/PQS/001.6.window.MergedG4.bed") %>% as.data.frame()
table(window$V1)
# chr2L chr2R chr3L chr3R  chr4  chrX  chrY 
# 2352  2529  2812  3208   135  2355   367
tapply(window$V4, window$V1, mean) 
# chr2L      chr2R      chr3L      chr3R       chr4       chrX       chrY 
# 1.18494898 1.30407276 1.16607397 1.33416459 0.11111111 1.18513800 0.01089918 
mean(window[window$V1!="chr4"&window$V1!="chrX"&window$V1!="chrY","V4"]) #G4在常染色体上平均密度是1.251628

has = fread("/home/yuss/flyG4/result/PQS/001.6.HAS.MergedG4.bed") %>% as.data.frame()
sum(has$V4)/nrow(has) # 0.8365759
nrow(has) #257

p = matrix(0,nrow = 5,ncol = 3)
for (condition in 1) {
  for (i in 1:5) {
    chr = paste0("chr",c("2L","2R","3L","3R","X"))[i]
    p[i,condition] = wilcox.test(window[window$V1==chr,condition+3],has[,condition+3])$p.value
  }
}
##has p.value=1.062922e-02

h4k = fread("/home/yuss/flyG4/result/PQS/001.6.H4K16ac.bothsex.merge.MergedG4s.bed") %>% as.data.frame()
h4k$len = h4k$V3 - h4k$V2 
sum(h4k$V4)/sum(h4k$len)*10^4 #1.233117
p = matrix(0,nrow = 5,ncol = 3)
for (condition in 1) {
  for (i in 1:5) {
    chr = paste0("chr",c("2L","2R","3L","3R","X"))[i]
    p[i,condition] = wilcox.test(window[window$V1==chr,condition+3],h4k[,condition+3])$p.value
  }
}
##h4k16ac p.value=8.491871e-299
df = data.frame(name = c("2L","2R","3L","3R","X","HAS","H4K16ac"),
                num = c(as.numeric(tapply(window$V4, window$V1, mean))[c(1:4,6)],sum(has$V4)/nrow(has),sum(h4k$V4)/sum(h4k$len)*10^4))
df$name %<>% factor(.,levels = c("2L","2R","3L","3R","X","HAS","H4K16ac"))
ggplot(df,aes(name,num))+geom_bar(stat = "identity",position = "dodge")+
  theme_bw()+ ylab("Density (No./kb)")+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=10,colour = "black",angle = 45,hjust = 1,vjust = 1),axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,colour = "black"),axis.title.y = element_text(size=12))+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(df$num)*1.1))

#*KC S2单独做-----------------------------------------------------------------------------------------
rm(list = ls());gc();rm(list = ls())
kcwindow = fread("/home/yuss/flyG4/result/PQS/001.6.window.kc.bed") %>% as.data.frame()
table(kcwindow$V1)
# chr2L chr2R chr3L chr3R  chr4  chrX  chrY 
# 2352  2529  2812  3208   135  2355   367
tapply(kcwindow$V4, kcwindow$V1, mean) 
# chr2L      chr2R      chr3L      chr3R       chr4       chrX       chrY 
# 0.23256803 0.24911032 0.22688478 0.31421446 0.00000000 0.67218684 0.01089918 
mean(kcwindow[kcwindow$V1!="chr4"&kcwindow$V1!="chrX"&kcwindow$V1!="chrY","V4"]) #G4在常染色体上平均密度是0.2589671

kc_has = fread("/home/yuss/flyG4/result/PQS/001.6.HAS.kc.bed") %>% as.data.frame()
sum(kc_has$V4)/nrow(kc_has) # 0.2840467
nrow(kc_has) #257

p = matrix(0,nrow = 5,ncol = 3)
for (condition in 1) {
  for (i in 1:5) {
    chr = paste0("chr",c("2L","2R","3L","3R","X"))[i]
    p[i,condition] = wilcox.test(kcwindow[kcwindow$V1==chr,condition+3],kc_has[,condition+3])$p.value
  }
}
##kc_has p.value=5.781405e-08

kc_h4k = fread("/home/yuss/flyG4/result/PQS/001.6.H4K16ac.bothsex.kc.bed") %>% as.data.frame()
kc_h4k$len = kc_h4k$V3 - kc_h4k$V2 
sum(kc_h4k$V4)/sum(kc_h4k$len)*10^4 #0.3949798
p = matrix(0,nrow = 5,ncol = 3)
for (condition in 1) {
  for (i in 1:5) {
    chr = paste0("chr",c("2L","2R","3L","3R","X"))[i]
    p[i,condition] = wilcox.test(kcwindow[kcwindow$V1==chr,condition+3],kc_h4k[,condition+3])$p.value
  }
}
##kc_h4k16ac p.value=0
df = data.frame(name = c("2L","2R","3L","3R","X","HAS","H4K16ac"),
                num = c(as.numeric(tapply(kcwindow$V4, kcwindow$V1, mean))[c(1:4,6)],sum(kc_has$V4)/nrow(kc_has),sum(kc_h4k$V4)/sum(kc_h4k$len)*10^4))
df$name %<>% factor(.,levels = c("2L","2R","3L","3R","X","HAS","H4K16ac"))
ggplot(df,aes(name,num))+geom_bar(stat = "identity",position = "dodge")+
  theme_bw()+ ylab("Density (No./kb)")+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=10,colour = "black",angle = 45,hjust = 1,vjust = 1),axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,colour = "black"),axis.title.y = element_text(size=12))+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(df$num)*1.1))

#*S2-----------------------------------------------------------------------------------------------------
s2window = fread("/home/yuss/flyG4/result/PQS/001.6.window.s2.bed") %>% as.data.frame()
table(s2window$V1)
# chr2L chr2R chr3L chr3R  chr4  chrX  chrY 
# 2352  2529  2812  3208   135  2355   367
tapply(s2window$V4, s2window$V1, mean) 
# chr2L      chr2R      chr3L      chr3R       chr4       chrX       chrY 
# 0.26530612 0.18149466 0.22261735 0.14869077 0.05925926 0.09893843 0.00000000 
mean(s2window[s2window$V1!="chr4"&s2window$V1!="chrX"&s2window$V1!="chrY","V4"]) #G4在常染色体上平均密度是0.2005321

s2_has = fread("/home/yuss/flyG4/result/PQS/001.6.HAS.s2.bed") %>% as.data.frame()
sum(s2_has$V4)/nrow(s2_has) # 0.1245136
nrow(s2_has) #257

p = matrix(0,nrow = 5,ncol = 3)
for (condition in 1) {
  for (i in 1:5) {
    chr = paste0("chr",c("2L","2R","3L","3R","X"))[i]
    p[i,condition] = wilcox.test(s2window[s2window$V1==chr,condition+3],s2_has[,condition+3])$p.value
  }
}
##s2_has p.value=0.1535648703

s2_h4k = fread("/home/yuss/flyG4/result/PQS/001.6.H4K16ac.bothsex.s2.bed") %>% as.data.frame()
s2_h4k$len = s2_h4k$V3 - s2_h4k$V2 
sum(s2_h4k$V4)/sum(s2_h4k$len)*10^4 #0.1579155
p = matrix(0,nrow = 5,ncol = 3)
for (condition in 1) {
  for (i in 1:5) {
    chr = paste0("chr",c("2L","2R","3L","3R","X"))[i]
    p[i,condition] = wilcox.test(s2window[s2window$V1==chr,condition+3],s2_h4k[,condition+3])$p.value
  }
}
##s2_h4k16ac p.value=1.068660e-16
df = data.frame(name = c("2L","2R","3L","3R","X","HAS","H4K16ac"),
                num = c(as.numeric(tapply(s2window$V4, s2window$V1, mean))[c(1:4,6)],sum(s2_has$V4)/nrow(s2_has),sum(s2_h4k$V4)/sum(s2_h4k$len)*10^4))
df$name %<>% factor(.,levels = c("2L","2R","3L","3R","X","HAS","H4K16ac"))
ggplot(df,aes(name,num))+geom_bar(stat = "identity",position = "dodge")+
  theme_bw()+ ylab("Density (No./kb)")+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=10,colour = "black",angle = 45,hjust = 1,vjust = 1),axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,colour = "black"),axis.title.y = element_text(size=12))+
  scale_y_continuous(expand = c(0,0),limits = c(0,max(df$num)*1.1))
