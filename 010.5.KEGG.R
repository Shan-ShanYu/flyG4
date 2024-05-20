rm(list = ls());gc();rm(list = ls())
Num = "010.5."
# devtools::install_github("YuLab-SMU/DOSE")
# devtools::install_github("YuLab-SMU/HDO.db")
# devtools::install_github('YuLab-SMU/clusterProfiler')  
library(clusterProfiler)
library(HDO.db)
library(DOSE)
library(dplyr)
library(data.table)
#### kc.Phen ####
kc.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.DE.kc.Phen.txt") %>% as.data.frame()
df <- kc.Phen[kc.Phen$group=="Up",]
kcPhenup <- unlist(df$geneid)

##富集分析
library(org.Dm.eg.db)
enrich_go_kegg<-function(gene_id){
  result<-list()
  ensembl_2_entrezid<-bitr(gene_id,fromType ="ENSEMBL", toType = c("ENTREZID","SYMBOL"), OrgDb = org.Dm.eg.db)
  ego_ALL <- enrichGO(gene = as.character(ensembl_2_entrezid$ENTREZID),
                      OrgDb=org.Dm.eg.db,
                      keyType = "ENTREZID",
                      ont = "ALL",
                      pAdjustMethod = "BH",
                      minGSSize = 1,
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05,
                      readable = TRUE)
  result$ego_all<-as.data.frame(ego_ALL)
  ensembl_2_kegg_id<-bitr_kegg(ensembl_2_entrezid$ENTREZID,fromType = "ncbi-geneid",toType = "kegg",organism = "dme")
  kegg<-enrichKEGG(gene = ensembl_2_entrezid$ENTREZID,organism = "dme", keyType = "ncbi-geneid",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
  result$kegg<-kegg@result
  result$ego_ALL<-ego_ALL
  return(result)
}

enrich_kcPhenup <- enrich_go_kegg(kcPhenup)
enrich_kcPhenup_ego <- enrich_kcPhenup$ego_all
enrich_kcPhenup_kegg <- enrich_kcPhenup$kegg
enrich_kcPhenup_ego_ALL <- enrich_kcPhenup$ego_ALL
dotplot(enrich_kcPhenup_ego_ALL,showCategory=10,split="ONTOLOGY",label_format=60,title="PhenDC3 up regulated(Kc) GO enrichment") + facet_grid(ONTOLOGY~., scale='free')
##showCategory指定展示的GO Terms的个数，默认展示显著富集的top10个
##label_format=60左边的名称一行显示60
##facet_grid(ONTOLOGY~., scale='free')将三个基因本体分开
library(ggplot2)
enrich_kcPhenup_kegg <- arrange(enrich_kcPhenup_kegg, enrich_kcPhenup_kegg$p.adjust)
enrich_kcPhenup_kegg$Description1 <- sapply(strsplit(enrich_kcPhenup_kegg$Description, " - "), function(x) x[1])
ggplot(enrich_kcPhenup_kegg[1:5,],aes(x=Count/87,y=Description1,colour=-1*log10(pvalue),size=Count))+
  geom_point()+
  scale_size(range = c(2, 8))+ #scale_size修改图中的点的大小，范围是2到8
  scale_color_gradient(low = "blue", high = "red")+
  theme_bw()+
  ylab("KEGG Pathway Terms")+
  xlab("Gene Ratio")+
  labs(color=expression(-log[10](PValue)), title = "PhenDC3 up regulated(Kc) kegg enrichment") +#expression函数改变样式，[]是用来添加下标，^是用来添加上标
  theme(text = element_text(size=15))

##down
df <- kc.Phen[kc.Phen$group=="Down",]
kcPhendown <- unlist(df$geneid)

##富集分析
enrich_kcPhendown <- enrich_go_kegg(kcPhendown)
enrich_kcPhendown_ego <- enrich_kcPhendown$ego_all
enrich_kcPhendown_kegg <- enrich_kcPhendown$kegg
enrich_kcPhendown_ego_ALL <- enrich_kcPhendown$ego_ALL
dotplot(enrich_kcPhendown_ego_ALL,showCategory=10,split="ONTOLOGY",label_format=60,title="PhenDC3 down regulated(Kc) GO enrichment") + facet_grid(ONTOLOGY~., scale='free')
##showCategory指定展示的GO Terms的个数，默认展示显著富集的top10个
##label_format=60左边的名称一行显示60
##facet_grid(ONTOLOGY~., scale='free')将三个基因本体分开
enrich_kcPhendown_kegg <- arrange(enrich_kcPhendown_kegg, enrich_kcPhendown_kegg$p.adjust)
enrich_kcPhendown_kegg$Description1 <- sapply(strsplit(enrich_kcPhendown_kegg$Description, " - "), function(x) x[1])
ggplot(enrich_kcPhendown_kegg[1:8,],aes(x=Count/64,y=Description1,colour=-1*log10(pvalue),size=Count))+
  geom_point()+
  scale_size(range = c(2, 8))+ #scale_size修改图中的点的大小，范围是2到8
  scale_color_gradient(low = "blue", high = "red")+
  theme_bw()+
  ylab("KEGG Pathway Terms")+
  xlab("Gene Ratio")+
  labs(color=expression(-log[10](PValue)), title = "PhenDC3 down regulated(Kc) kegg enrichment") +#expression函数改变样式，[]是用来添加下标，^是用来添加上标
  theme(text = element_text(size=15))

rm(list = ls());gc();rm(list = ls())
library(clusterProfiler)
library(HDO.db)
library(DOSE)
library(dplyr)
library(data.table)
#### s2.Phen ####
s2.Phen <- fread("/home/yuss/flyG4/result/KcS2.RNAseq/010.1.DE.s2.Phen.txt") %>% as.data.frame()
df <- s2.Phen[s2.Phen$group=="Up",]
s2Phenup <- unlist(df$geneid)

##富集分析
library(org.Dm.eg.db)
enrich_go_kegg<-function(gene_id){
  result<-list()
  ensembl_2_entrezid<-bitr(gene_id,fromType ="ENSEMBL", toType = c("ENTREZID","SYMBOL"), OrgDb = org.Dm.eg.db)
  ego_ALL <- enrichGO(gene = as.character(ensembl_2_entrezid$ENTREZID),
                      OrgDb=org.Dm.eg.db,
                      keyType = "ENTREZID",
                      ont = "ALL",
                      pAdjustMethod = "BH",
                      minGSSize = 1,
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05,
                      readable = TRUE)
  result$ego_all<-as.data.frame(ego_ALL)
  ensembl_2_kegg_id<-bitr_kegg(ensembl_2_entrezid$ENTREZID,fromType = "ncbi-geneid",toType = "kegg",organism = "dme")
  kegg<-enrichKEGG(gene = ensembl_2_entrezid$ENTREZID,organism = "dme", keyType = "ncbi-geneid",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
  result$kegg<-kegg@result
  result$ego_ALL<-ego_ALL
  return(result)
}

enrich_s2Phenup <- enrich_go_kegg(s2Phenup)
enrich_s2Phenup_ego <- enrich_s2Phenup$ego_all
enrich_s2Phenup_kegg <- enrich_s2Phenup$kegg
enrich_s2Phenup_ego_ALL <- enrich_s2Phenup$ego_ALL
dotplot(enrich_s2Phenup_ego_ALL,showCategory=10,split="ONTOLOGY",label_format=60,title="PhenDC3 up regulated(S2) GO enrichment") + facet_grid(ONTOLOGY~., scale='free')

library(ggplot2)
enrich_s2Phenup_kegg <- arrange(enrich_s2Phenup_kegg, enrich_s2Phenup_kegg$p.adjust)
enrich_s2Phenup_kegg$Description1 <- sapply(strsplit(enrich_s2Phenup_kegg$Description, " - "), function(x) x[1])
ggplot(enrich_s2Phenup_kegg[1:11,],aes(x=Count/96,y=Description1,colour=-1*log10(pvalue),size=Count))+
  geom_point()+
  scale_size(range = c(2, 8))+ #scale_size修改图中的点的大小，范围是2到8
  scale_color_gradient(low = "blue", high = "red")+
  theme_bw()+
  ylab("KEGG Pathway Terms")+
  xlab("Gene Ratio")+
  labs(color=expression(-log[10](PValue)), title = "PhenDC3 up regulated(S2) kegg enrichment") +#expression函数改变样式，[]是用来添加下标，^是用来添加上标
  theme(text = element_text(size=15))

##down
df <- s2.Phen[s2.Phen$group=="Down",]
s2Phendown <- unlist(df$geneid)

##富集分析
enrich_s2Phendown <- enrich_go_kegg(s2Phendown)
enrich_s2Phendown_ego <- enrich_s2Phendown$ego_all
enrich_s2Phendown_kegg <- enrich_s2Phendown$kegg
enrich_s2Phendown_ego_ALL <- enrich_s2Phendown$ego_ALL
dotplot(enrich_s2Phendown_ego_ALL,showCategory=10,split="ONTOLOGY",label_format=60,title="PhenDC3 down regulated(S2) GO enrichment") + facet_grid(ONTOLOGY~., scale='free')
##showCategory指定展示的GO Terms的个数，默认展示显著富集的top10个
##label_format=60左边的名称一行显示60
##facet_grid(ONTOLOGY~., scale='free')将三个基因本体分开
enrich_s2Phendown_kegg <- arrange(enrich_s2Phendown_kegg, enrich_s2Phendown_kegg$p.adjust)
enrich_s2Phendown_kegg$Description1 <- sapply(strsplit(enrich_s2Phendown_kegg$Description, " - "), function(x) x[1])
ggplot(enrich_s2Phendown_kegg[1:15,],aes(x=Count/94,y=Description1,colour=-1*log10(pvalue),size=Count))+
  geom_point()+
  scale_size(range = c(2, 8))+ #scale_size修改图中的点的大小，范围是2到8
  scale_color_gradient(low = "blue", high = "red")+
  theme_bw()+
  ylab("KEGG Pathway Terms")+
  xlab("Gene Ratio")+
  labs(color=expression(-log[10](PValue)), title = "PhenDC3 down regulated(S2) kegg enrichment") +#expression函数改变样式，[]是用来添加下标，^是用来添加上标
  theme(text = element_text(size=15))
