rm(list = ls());gc();rm(list = ls())#清空
Num = "003.6."

# 4082 003.6.kc.2-complex.fasta
# 1652 003.6.kc.long-loop.fasta
# 98 003.6.kc.loop1-3.fasta
# 210 003.6.kc.loop4-5.fasta
# 340 003.6.kc.loop6-7.fasta
# 10172 003.6.kc.other.fasta
# 11408 003.6.kc.simple-bulge.fasta
# 27962 total

# 重新分类后的数目
# 4082 003.6.kc.2-complex.fasta
# 86 003.6.kc.long-loop.fasta
# 2444 003.6.kc.loop1-3.fasta
# 742 003.6.kc.loop4-5.fasta
# 296 003.6.kc.loop6-7.fasta
# 10170 003.6.kc.other.fasta
# 10142 003.6.kc.simple-bulge.fasta
# 27962 total

kc.data <- data.frame(
  group = c("loop1-3", "loop4-5", "loop6-7", "long-loop", "simple-bulge", "2-tetrad/complex bulges", "other"),
  num = c(2444/2, 742/2, 296/2, 86/2, 10142/2, 4082/2, 10170/2)
)
kc.data$percent <- kc.data$num/13981
kc.data$group <- factor(kc.data$group, levels = c("other","2-tetrad/complex bulges","simple-bulge","long-loop","loop6-7","loop4-5","loop1-3"))
colors <- c("#A3B7D9","#E4964F", "#4FAAB6", "#856FA0", "#9BB25F", "#BA5B54","#5786B6")
library(ggplot2)
ggplot(kc.data, aes(x="", y=num, fill=group))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y") + #原始饼图
  scale_fill_manual(values = colors, guide = guide_legend(reverse = T)) +
  geom_text(aes(label = paste0(round(percent * 100, 2), "%")), 
            position = position_stack(vjust = 0.5)) +
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank()##删除图例标题名称
  )

# s2
# 3434 003.6.s2.2-complex.fasta
# 78 003.6.s2.long-loop.fasta
# 2162 003.6.s2.loop1-3.fasta
# 640 003.6.s2.loop4-5.fasta
# 252 003.6.s2.loop6-7.fasta
# 8828 003.6.s2.other.fasta
# 8614 003.6.s2.simple-bulge.fasta
# 24008 total
s2.data <- data.frame(
  group = c("loop1-3", "loop4-5", "loop6-7", "long-loop", "simple-bulge", "2-tetrad/complex bulges", "other"),
  num = c(2162/2, 640/2, 252/2, 78/2, 8614/2, 3434/2, 8828/2)
)
s2.data$percent <- s2.data$num/12004
s2.data$group <- factor(s2.data$group, levels = c("other","2-tetrad/complex bulges","simple-bulge","long-loop","loop6-7","loop4-5","loop1-3"))
colors <- c("#A3B7D9","#E4964F", "#4FAAB6", "#856FA0", "#9BB25F", "#BA5B54","#5786B6")
ggplot(s2.data, aes(x="", y=num, fill=group))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y") + #原始饼图
  scale_fill_manual(values = colors, guide = guide_legend(reverse = T)) +
  geom_text(aes(label = paste0(round(percent * 100, 2), "%")), 
            position = position_stack(vjust = 0.5)) +
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank()##删除图例标题名称
  )

