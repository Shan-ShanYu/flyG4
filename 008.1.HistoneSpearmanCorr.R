rm(list = ls());gc();rm(list = ls())
Num = "008.1."
##install.packages("ggcorrplot")
library(ggcorrplot)
s2.SpearmanCorr <- fread('/home/yuss/flyG4/result/Histone.modifications/s2.scatterplot.SpearmanCorr.readCounts.tab') %>% data.frame()
row.names(s2.SpearmanCorr) <- c("S2-G4-CUTTag","ATAC-seq","POLR2A", "H3K27ac","H3K4me1","H3K4me3","H4k16ac")
s2.SpearmanCorr <- s2.SpearmanCorr[,-1]
colnames(s2.SpearmanCorr) <- c("S2-G4-CUTTag","ATAC-seq","POLR2A", "H3K27ac","H3K4me1","H3K4me3","H4k16ac")
s2.SpearmanCorr <- as.data.frame(sapply(s2.SpearmanCorr, as.numeric))
row.names(s2.SpearmanCorr) <- c("S2-G4-CUTTag","ATAC-seq","POLR2A", "H3K27ac","H3K4me1","H3K4me3","H4k16ac")
ggcorrplot(s2.SpearmanCorr, method = "circle", type = "upper", ggtheme = ggplot2::theme_classic(), title = "", 
           show.legend = TRUE, legend.title = "Corr", show.diag = T, 
           colors = c("blue", "white", "red"), outline.color = "grey", 
           # hc.order = T, hc.method = "complete",
           lab = FALSE, lab_col = "black", 
           lab_size = 4, p.mat = NULL, sig.level = 0.05, insig = c("pch", "blank"), pch = 4, pch.col = "black", pch.cex = 5, tl.cex = 12, 
           tl.col = "black", tl.srt = 45, digits = 2)

##另外一种画法
s2.SpearmanCorr <- fread('/home/yuss/flyG4/result/Histone.modifications/s2.scatterplot.SpearmanCorr.readCounts.tab')
row.names(s2.SpearmanCorr) <- c("S2-G4-CUTTag","ATAC-seq","POLR2A", "H3K27ac","H3K4me1","H3K4me3","H4k16ac")
s2.SpearmanCorr <- s2.SpearmanCorr[,-1]
colnames(s2.SpearmanCorr) <- c("S2-G4-CUTTag","ATAC-seq","POLR2A", "H3K27ac","H3K4me1","H3K4me3","H4k16ac")
class(s2.SpearmanCorr)
# s2.SpearmanCorr <- as.numeric(as.matrix(s2.SpearmanCorr))
s2.SpearmanCorr <- as.matrix(sapply(s2.SpearmanCorr, as.numeric))
row.names(s2.SpearmanCorr) <- c("S2-G4-CUTTag","ATAC-seq","POLR2A", "H3K27ac","H3K4me1","H3K4me3","H4k16ac")
corrplot(s2.SpearmanCorr)
pdf("/home/yuss/flyG4/result/Histone.modifications/Picture/008.1.S2HistoneSpearmanCorr.pdf",width = 5,height = 5)
corrplot(s2.SpearmanCorr, method = "circle", type = "upper",
         tl.col = "black", tl.cex = 0.8, tl.srt = 45,tl.pos = "lt")
corrplot(s2.SpearmanCorr, method = "number", type = "lower",
         tl.col = "n", tl.cex = 0.8, tl.pos = "n",
         add = T)
dev.off()
##ggsave(filename = paste0("/home/yuss/flyG4/result/Histone.modifications/Picture/",Num,"S2HistoneSpearmanCorr.pdf"),
device = "pdf",width = 4,height = 4)

kc.SpearmanCorr <- fread('/home/yuss/flyG4/result/Histone.modifications/kc.scatterplot.SpearmanCorr.readCounts.tab')
row.names(kc.SpearmanCorr) <- c("Kc-G4-CUTTag","ATAC-seq","POLR2A", "H3K27ac","H3K4me1","H3K4me3","H4k16ac")
kc.SpearmanCorr <- kc.SpearmanCorr[,-1]
colnames(kc.SpearmanCorr) <- c("Kc-G4-CUTTag","ATAC-seq","POLR2A", "H3K27ac","H3K4me1","H3K4me3","H4k16ac")
class(kc.SpearmanCorr)
kc.SpearmanCorr <- as.matrix(sapply(kc.SpearmanCorr, as.numeric))
row.names(kc.SpearmanCorr) <- c("Kc-G4-CUTTag","ATAC-seq","POLR2A", "H3K27ac","H3K4me1","H3K4me3","H4k16ac")
corrplot(kc.SpearmanCorr)
pdf("/home/yuss/flyG4/result/Histone.modifications/Picture/008.1.KcHistoneSpearmanCorr.pdf",width = 5,height = 5)
corrplot(kc.SpearmanCorr, method = "circle", type = "upper",
         tl.col = "black", tl.cex = 0.8, tl.srt = 45,tl.pos = "lt")
corrplot(kc.SpearmanCorr, method = "number", type = "lower",
         tl.col = "n", tl.cex = 0.8, tl.pos = "n",
         add = T)
dev.off()
##ggsave(filename = paste0("/home/yuss/flyG4/result/Histone.modifications/Picture/",Num,"S2HistoneSpearmanCorr.pdf"),
device = "pdf",width = 4,height = 4)