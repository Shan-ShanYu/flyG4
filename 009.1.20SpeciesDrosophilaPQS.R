#生成20种果蝇的pqs(分数)
#### 生成score>50的bed ####
rm(list = ls());gc();rm(list = ls())
Num = "009.1."

library(tidyverse)
library(pqsfinder)
library(Biostrings)
#a输入基因组序列fa文件路径，b输出pqsbed文件路径
a <- "/home/yuss/flyG4/data/Drosophila_genome/Dmel.fna"
b <- "/home/yuss/flyG4/result/20SpeciesDrosophilaPQS/009.1.Dmel.pqs.bed"
#将基因组序列文件导入R成为变量genome
genome <- readDNAStringSet(a)
#使用pqsfinder在每个染色体上分别注释pqs，结果赋给chr_pqs。再从得到的结果中提取所需数据成为数据框chr?_df，包括：各个PQS的start，end，strand，score
allpqs <- data.frame()
for (i in c(1:length(genome@ranges@NAMES))) {
  chri_pqs <- pqsfinder(genome[[i]],overlapping = FALSE,min_score = 50)
  chri_df <- list(chri_pqs@ranges,chri_pqs@elementMetadata)%>%
    data.frame()%>%
    mutate(chr=paste("chr",strsplit(genome@ranges@NAMES[i]," ")[[1]][1],sep = ""))%>%
    select(15,1,2,5,4)
  allpqs <- rbind(allpqs,chri_df)
}
#产生id
allpqs <- mutate(allpqs,PQSID=paste("id",n=c(1:nrow(allpqs)),sep="_"))%>%
  select(1,2,3,6,4,5)
#储存
data.table::fwrite(allpqs,file =b,sep = '\t',row.names = F,quote = F,col.names = F)

#### 生成score>55的bed ####
rm(list = ls());gc();rm(list = ls())
Num = "009.1."
setwd("/home/yuss/flyG4/result/20SpeciesDrosophilaPQS")
files <- list.files("/home/yuss/flyG4/result/20SpeciesDrosophilaPQS",pattern = "\\.bed$")
path <- '/home/yuss/flyG4/result/20SpeciesDrosophilaPQS'
filepath <- sapply(files, function(x){
  paste(path,x,sep = '/')
})
data <- list()
for (i in 1:length(files)){
  data[[i]] <- fread(filepath[[i]])
}
##分开的文件
# for (i in 1:length(files)){
#   var_name <- gsub('pqs.bed', '',files[i])
#   var_name <- gsub('009.1.', '',var_name)
#   assign(var_name, read.table(files[i], sep = '\t', header = F)) ##assign()函数将一个读取的数据框对象分配给先前定义的变量名 var_name
# }

# 创建一个新的存储路径
output_path <- '/home/yuss/flyG4/result/20SpeciesDrosophilaPQS'
dir.create(output_path, showWarnings = FALSE)

# 循环筛选和保存
for (i in seq_along(data)) {
  # 筛选第五列大于55的数据
  filtered_data <- data[[i]] %>% filter(V5 > 55)
  
  # 构造新的文件名
  new_filename <- gsub('.bed', '.55bed',files[i])
  new_filename <- paste0(output_path, '/', new_filename)
  
  # 保存筛选后的数据
  write.table(filtered_data, file = new_filename, sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
}

#### 生成score>60的bed ####
rm(list = ls());gc();rm(list = ls())
Num = "009.1."
setwd("/home/yuss/flyG4/result/20SpeciesDrosophilaPQS")
files <- list.files("/home/yuss/flyG4/result/20SpeciesDrosophilaPQS",pattern = "\\.bed$")
path <- '/home/yuss/flyG4/result/20SpeciesDrosophilaPQS'
filepath <- sapply(files, function(x){
  paste(path,x,sep = '/')
})
data <- list()
for (i in 1:length(files)){
  data[[i]] <- fread(filepath[[i]])
}

# 创建一个新的存储路径
output_path <- '/home/yuss/flyG4/result/20SpeciesDrosophilaPQS'
dir.create(output_path, showWarnings = FALSE)

# 循环筛选和保存
for (i in seq_along(data)) {
  # 筛选第五列大于55的数据
  filtered_data <- data[[i]] %>% filter(V5 > 60)
  
  # 构造新的文件名
  new_filename <- gsub('.bed', '.60bed',files[i])
  new_filename <- paste0(output_path, '/', new_filename)
  
  # 保存筛选后的数据
  write.table(filtered_data, file = new_filename, sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
}

