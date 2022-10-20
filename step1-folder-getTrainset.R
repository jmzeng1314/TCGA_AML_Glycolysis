rm(list=ls())
options(stringsAsFactors = F)
library(stringr)  
library(data.table) 
library(tidyverse)

# 0. create folder of task----
dir.create('output',recursive = T)
dir.create('output/plot',recursive = T)
dir.create('output/rdata',recursive = T)
dir.create('output/text',recursive = T)

# 1. tidy clinical information----
phenotype.org <- fread('input/TCGA-LAML.GDC_phenotype.tsv.gz',
                      data.table = F )
phenotype.org[1:6, 1:6]
dim(phenotype.org)
colnames(phenotype.org)[41]
table(phenotype.org$leukemia_french_american_british_morphology_code)

## 1.1 remove the phenotype information of M3
phenotype <- phenotype.org %>% 
  filter(leukemia_french_american_british_morphology_code != 'M3')
phenotype[1:4, 1:4]
dim(phenotype)


# 2. tidy survival data----
surv.org <- fread('input/TCGA-LAML.survival.tsv',data.table = F)
head(surv.org)
dim(surv.org)
colnames(surv.org)
## filter time
surv <- surv.org %>%
  filter(OS.time >= 30) %>%
  select(sample, OS, OS.time)
head(surv)
dim(surv)
## days to year
surv$time <- surv$OS.time/365
hist(surv$time)
surv <- surv[,-3]
head(surv)
colnames(surv)[2] <- 'event'
head(surv)


# 3. tidy expression ----
## dataset: gene expression RNAseq - HTSeq - Counts 
## unit: log2(count+1)
a <- fread('input/TCGA-LAML.htseq_counts.tsv.gz',data.table = F) 
head(a[ ,1:4])
tail(a[ ,1:4])

mat <- a[1:60483,]
mat = a
rownames(mat) <- a$Ensembl_ID  # keep point
mat[1:4,1:4]
mat <- mat[,-1]
mat[1:4,1:4]


## 2.1 tidy data to integer----
org_dat <- floor(2^mat - 1 ) 
org_dat[1:4,1:4] 

## 2.2 filter data----
n <- floor(ncol(org_dat)/20)
n
keep_feature <- rowSums (org_dat > 2) > n
table(keep_feature)

mat <- org_dat[keep_feature, ]
mat <- mat[, colSums(mat) > 1000000]  # !!note: can change to get more genes!!
mat[1:4,1:4]
dim(mat) 
head(colnames(mat))

## 2.3 ensg to symbol ----
library(AnnoProbe)
b <- fread('input/gencode.v22.annotation.gene.probeMap',header = T, data.table = F)
head(b)  # ensg has point
# 放弃这个方法转ID，因为失败率太高。
tmp= annoGene(b$gene,'SYMBOL')
tail(sort(table(tmp$biotypes)))


d <- fread('input/ensembID2type.txt', header = F, data.table = F)
head(d)
b <- merge(b,d,by.x='id',by.y='V1')
head(b)
length(unique(b[match(rownames(mat),b$id),2]))

## 2.4 get unique symbol----
dat <- mat
ids <-  b[match(rownames(dat),b$id),1:2]
head(ids)

colnames(ids) <- c('probe_id','symbol')
ids <- ids[ids$symbol != '',]
ids <- ids[ids$probe_id %in%  rownames(dat),]
dat[1:4,1:4]
dat <- dat[ids$probe_id,]

ids$median <- apply(dat,1,median) #ids新建median这一列，列名为median，同时对dat这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
ids <- ids[order(ids$symbol,ids$median,decreasing = T),]#对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
ids <- ids[!duplicated(ids$symbol),]#将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果s
dat <- dat[ids$probe_id,] #新的ids取出probe_id这一列，将dat按照取出的这一列中的每一行组成一个新的dat
rownames(dat) <- ids$symbol#把ids的symbol这一列中的每一行给dat作为dat的行名
dat[1:4,1:4]  #保留每个基因ID第一次出现的信息

## 2.5 check housekeeping gene----
mat=dat
mat[1:4,1:4]
fivenum(mat['GAPDH',])
fivenum(mat['ACTB',])

## 2.6 check sample id---- 
tmp = as.numeric( substring(colnames(mat),14,15)) < 10
table(tmp)
nchar(head(colnames(mat)))
#colnames(mat) = substring(colnames(mat),1,16)

## 2.7 divide into protein coding and non protein coding----
tp=b[match(rownames(mat),b$gene),7]
tail(sort(table(tp)))

pd_mat=mat[tp=='protein_coding',]
pd_mat[1:4, 1:4]
dim(pd_mat)

## special filter 
exp = pd_mat[apply(pd_mat, 1, function(x) sum(x > 0) > 0.5*ncol(pd_mat)), ]
exp[1:4, 1:4]
dim(exp)

# 4. overlap sample id----
kp <- Reduce(intersect, list(colnames(exp), 
                                 phenotype$submitter_id.samples,
                                 surv$sample))
head(kp)

mat <- exp[,kp]
mat[1:4, 1:4]
pheno_data <- phenotype[match(colnames(mat), phenotype$submitter_id.samples),]
head(pheno_data)[,1:4]
dim(pheno_data)
surv_data <- surv[match(colnames(mat), surv$sample),]
head(surv_data)
## check id order 
identical(colnames(mat), pheno_data$submitter_id.samples)
identical(colnames(mat), surv_data$sample)

# 5. target genes----
## glycolysis
genes_target <- read.table('input/Glycolysis_gene_list.txt')[,1]
head(genes_target)
length(unique(genes_target))

# 6. save data----
save(mat, pheno_data,surv_data,genes_target, 
     file = 'output/rdata/prepare_data_for_cox.Rdata')




