rm(list=ls())
options(stringsAsFactors = F)

library(tidyverse)
library(tidyr)
library(stringr)
library(data.table)

# 1. load data----
library(AnnoProbe)
library(GEOquery) 
## test dataset basic information
testset <- data.frame(dataset = paste0('dateset',1:5),
                      GSE = c(rep('GSE37642',2),rep('GSE12417',2), 'GSE71014'),
                      GPL = c(rep(c('GPL570','GPL96'),2),'GPL10558'))
head(testset)


## 1.1 test dataset 5----
gse_number <- testset$GSE[5]
gse_number

# gset <- geoChina(gse_number)
load('input/GSE71014_eSet.Rdata')

gset
gset[[1]]
a=gset[[1]] 
dat=exprs(a)
dat[1:4, 1:4]
dim(dat)
boxplot(dat[,1:4],las=2)

# dat=limma::normalizeBetweenArrays(dat)
# boxplot(dat[,1:4],las=2)
# ids <- idmap('GPL10558','soft')
load('input/GPL10558_soft.rda')

head(ids)
ids=ids[ids$symbol != '',]
dat=dat[rownames(dat) %in% ids$ID,]
ids=ids[match(rownames(dat),ids$ID),]
head(ids) 
colnames(ids)=c('probe_id','symbol')  
ids$probe_id=as.character(ids$probe_id)
head(ids)
rownames(dat)=ids$probe_id
dat[1:4,1:4] 

## gene symbols
ids=ids[ids$probe_id %in%  rownames(dat),]
dat[1:4,1:4]   
dat=dat[ids$probe_id,] 
ids$median=apply(dat,1,median)
ids=ids[order(ids$symbol,ids$median,decreasing = T),]
ids=ids[!duplicated(ids$symbol),]
dat=dat[ids$probe_id,]
rownames(dat)=ids$symbol
dat[1:4,1:4]
fivenum(dat['ACTB',])
fivenum(dat['GAPDH',])

## survival data
pd=pData(a) 
colnames(pd) ## no gender, no age
head(pd)[,1:4]
surv_data <- pd[,c(35,36)]
head(surv_data)
colnames(surv_data) <- c('event', 'time')
surv_data$event <- as.integer(str_sub(surv_data$event, 17))
table(surv_data$event)
str(surv_data)

## month to year
surv_data$time <- as.numeric(surv_data$time)/12
surv_data <- surv_data %>%
  filter(time >= 30/365)
hist(surv_data$time)
head(surv_data)
str(surv_data)

## sample id
identical(colnames(dat), rownames(surv_data))
kp <- intersect(colnames(dat), rownames(surv_data))
length(kp)
surv_data <- surv_data[match(kp,rownames(surv_data)),]
head(surv_data)
dat <- dat[,match(kp,colnames(dat))]
dat[1:4, 1:4]
identical(colnames(dat), rownames(surv_data))
testset5 <- list(expression = dat,
                 survival = surv_data)


## 1.2 test dataset 1&2----
gse_number <- testset$GSE[1]
gse_number
# gset <- geoChina(gse_number)
load('input/GSE37642_eSet.Rdata')

length(gset)

### 1.2.1 gpl570----
gset[[1]]
a=gset[[1]] 
dat=exprs(a)
dat[1:4, 1:4]
dim(dat)
boxplot(dat[,1:4],las=2)
# dat=limma::normalizeBetweenArrays(dat)
# boxplot(dat[,1:4],las=2)
# ids <- idmap('GPL570','soft')
load('input/GPL570_soft.rda')
head(ids)
ids=ids[ids$symbol != '',]
dat=dat[rownames(dat) %in% ids$ID,]
ids=ids[match(rownames(dat),ids$ID),]
head(ids) 
colnames(ids)=c('probe_id','symbol')  
ids$probe_id=as.character(ids$probe_id)
head(ids)
rownames(dat)=ids$probe_id
dat[1:4,1:4] 
## gene symbols
ids=ids[ids$probe_id %in%  rownames(dat),]
dat[1:4,1:4]   
dat=dat[ids$probe_id,] 
ids$median=apply(dat,1,median)
ids=ids[order(ids$symbol,ids$median,decreasing = T),]
ids=ids[!duplicated(ids$symbol),]
dat=dat[ids$probe_id,]
rownames(dat)=ids$symbol
dat[1:4,1:4]
dim(dat)
fivenum(dat['ACTB',])
fivenum(dat['GAPDH',])
## survival data
pd=pData(a) 
dim(pd)
colnames(pd) ## no gender
table(pd$`fab:ch1`)
head(pd)[,1:4]
surv_data <- pd %>%
  filter(`fab:ch1` != '3' & `fab:ch1` !='3v' & `fab:ch1` != 'NA')  %>%
  select(`age:ch1`, `life status:ch1`,`overall survival (days):ch1`)
dim(surv_data)
head(surv_data)
colnames(surv_data) <- c('age','status','time')
hist(as.numeric(surv_data$time))
table(surv_data$status)
surv_data$event <- as.integer(ifelse(surv_data$status == 'alive', 0, 1))
table(surv_data$event)
## day to year
surv_data$time <- as.numeric(surv_data$time)/365
surv_data <- surv_data %>%
  filter(time >= 30/356)
hist(surv_data$time)
head(surv_data)
str(surv_data)
## feature
surv_data$age <- as.numeric(surv_data$age)
head(surv_data)
str(surv_data)
## sample id
identical(colnames(dat), rownames(surv_data))
kp <- intersect(colnames(dat), rownames(surv_data))
length(kp)
surv_data <- surv_data[match(kp,rownames(surv_data)),]
head(surv_data)
dat <- dat[,match(kp,colnames(dat))]
dat[1:4, 1:4]
identical(colnames(dat), rownames(surv_data))

testset1 <- list(expression = dat,
                 survival = surv_data[,c('event','time')],
                 feature = surv_data[,c('event','time','age')])

### 1.2.2 gpl96----
gset[[2]]
a=gset[[2]] 
dat=exprs(a)
dat[1:4, 1:4]
dim(dat)
boxplot(dat[,1:4],las=2)
# dat=limma::normalizeBetweenArrays(dat)
# boxplot(dat[,1:4],las=2)
# ids <- idmap('GPL96','soft')
load('input/GPL96_soft.rda')
head(ids)
ids=ids[ids$symbol != '',]
dat=dat[rownames(dat) %in% ids$ID,]
ids=ids[match(rownames(dat),ids$ID),]
head(ids) 
colnames(ids)=c('probe_id','symbol')  
ids$probe_id=as.character(ids$probe_id)
head(ids)
rownames(dat)=ids$probe_id
dat[1:4,1:4] 
## gene symbols
ids=ids[ids$probe_id %in%  rownames(dat),]
dat[1:4,1:4]   
dat=dat[ids$probe_id,] 
ids$median=apply(dat,1,median)
ids=ids[order(ids$symbol,ids$median,decreasing = T),]
ids=ids[!duplicated(ids$symbol),]
dat=dat[ids$probe_id,]
rownames(dat)=ids$symbol
dat[1:4,1:4]
dim(dat)
fivenum(dat['ACTB',])
fivenum(dat['GAPDH',])
## survival data
pd=pData(a) 
dim(pd)
colnames(pd) ## no gender
table(pd$`fab:ch1`)
head(pd)[,1:4]
surv_data <- pd %>%
  filter(`fab:ch1` != '3' & `fab:ch1` !='3v' & `fab:ch1` != 'NA')  %>%
  select(`age:ch1`, `life status:ch1`,`overall survival (days):ch1`)
dim(surv_data)
head(surv_data)
colnames(surv_data) <- c('age','status','time')
hist(as.numeric(surv_data$time))
table(surv_data$status)
surv_data$event <- as.integer(ifelse(surv_data$status == 'alive', 0, 1))
table(surv_data$event)
## day to year
surv_data$time <- as.numeric(surv_data$time)/365
surv_data <- surv_data %>%
  filter(time >= 30/356)
hist(surv_data$time)
head(surv_data)
str(surv_data)
## feature
surv_data$age <- as.numeric(surv_data$age)
head(surv_data)
str(surv_data)
## sample id
identical(colnames(dat), rownames(surv_data))
kp <- intersect(colnames(dat), rownames(surv_data))
length(kp)
surv_data <- surv_data[match(kp,rownames(surv_data)),]
head(surv_data)
dat <- dat[,match(kp,colnames(dat))]
dat[1:4, 1:4]
identical(colnames(dat), rownames(surv_data))
testset2 <- list(expression = dat,
                 survival = surv_data[,c('event','time')],
                 feature = surv_data[,c('event','time','age')])

## 1.3 test dataset 3&4----
gse_number <- testset$GSE[3]
gse_number
# gset <- geoChina(gse_number)
load('input/GSE12417_eSet.Rdata')
gset
### 1.3.1 gpl570----
gset[[1]]
a=gset[[1]] 
dat=exprs(a)
dat[1:4, 1:4]
dim(dat)
boxplot(dat[,1:4],las=2)
# dat=limma::normalizeBetweenArrays(dat)
# boxplot(dat[,1:4],las=2)
# ids <- idmap('GPL570','soft')
load('input/GPL570_soft.rda')
head(ids)
ids=ids[ids$symbol != '',]
dat=dat[rownames(dat) %in% ids$ID,]
ids=ids[match(rownames(dat),ids$ID),]
head(ids) 
colnames(ids)=c('probe_id','symbol')  
ids$probe_id=as.character(ids$probe_id)
head(ids)
rownames(dat)=ids$probe_id
dat[1:4,1:4] 
## gene symbols
ids=ids[ids$probe_id %in%  rownames(dat),]
dat[1:4,1:4]   
dat=dat[ids$probe_id,] 
ids$median=apply(dat,1,median)
ids=ids[order(ids$symbol,ids$median,decreasing = T),]
ids=ids[!duplicated(ids$symbol),]
dat=dat[ids$probe_id,]
rownames(dat)=ids$symbol
dat[1:4,1:4]
dim(dat)
fivenum(dat['ACTB',])
fivenum(dat['GAPDH',])
## survival data
pd=pData(a) 
dim(pd)
colnames(pd)[1:33] ## no gender
head(pd)[,1:4]
surv_all <- str_split(pd$characteristics_ch1, pattern = ';',simplify = T) %>% as.data.frame()
head(surv_all)
dim(surv_all)
colnames(surv_all) <- c('fab','age','time','status')
rownames(surv_all) <- rownames(pd)
head(surv_all)
table(surv_all$fab) ## no M3

surv_aml <- surv_all %>%
  filter(fab != 'MDS normal karyotype (test set) MDS RAEB')
head(surv_aml)
surv_data <- data.frame(age = readr::parse_number(surv_aml$age),
                        time = readr::parse_number(surv_aml$time),
                        event = as.integer(str_split(surv_aml$status, ":", simplify = T)[,2]))
head(surv_data)
str(surv_data)
rownames(surv_data) <- rownames(surv_aml)
hist(as.numeric(surv_data$time))
table(surv_data$event)
## day to year
surv_data <- surv_data %>%
  filter(time >= 30)
surv_data$time <- as.numeric(surv_data$time)/365
hist(surv_data$time)
head(surv_data)
str(surv_data)
## sample id
identical(colnames(dat), rownames(surv_data))
kp <- intersect(colnames(dat), rownames(surv_data))
length(kp)
surv_data <- surv_data[match(kp,rownames(surv_data)),]
head(surv_data)
dat <- dat[,match(kp,colnames(dat))]
dat[1:4, 1:4]
identical(colnames(dat), rownames(surv_data))

testset3 <- list(expression = dat,
                 survival = surv_data[,c('event','time')],
                 feature = surv_data[,c('event','time','age')])
### 1.3.2 gpl96----
gset[[2]]
a=gset[[2]] 
dat=exprs(a)
dat[1:4, 1:4]
dim(dat)
boxplot(dat[,1:4],las=2)
# dat=limma::normalizeBetweenArrays(dat)
# boxplot(dat[,1:4],las=2)
# ids <- idmap('GPL96','soft')
load('input/GPL96_soft.rda')
head(ids)
ids=ids[ids$symbol != '',]
dat=dat[rownames(dat) %in% ids$ID,]
ids=ids[match(rownames(dat),ids$ID),]
head(ids) 
colnames(ids)=c('probe_id','symbol')  
ids$probe_id=as.character(ids$probe_id)
head(ids)
rownames(dat)=ids$probe_id
dat[1:4,1:4] 
## gene symbols
ids=ids[ids$probe_id %in%  rownames(dat),]
dat[1:4,1:4]   
dat=dat[ids$probe_id,] 
ids$median=apply(dat,1,median)
ids=ids[order(ids$symbol,ids$median,decreasing = T),]
ids=ids[!duplicated(ids$symbol),]
dat=dat[ids$probe_id,]
rownames(dat)=ids$symbol
dat[1:4,1:4]
dim(dat)
fivenum(dat['ACTB',])
fivenum(dat['GAPDH',])
## survival data
pd=pData(a) 
dim(pd)
colnames(pd)[1:4] ## no gender
head(pd)[,1:4]
surv_all <- str_split(pd$characteristics_ch1, pattern = ';',simplify = T) %>% as.data.frame()
head(surv_all)
dim(surv_all)
colnames(surv_all) <- c('fab','age','time','status')
rownames(surv_all) <- rownames(pd)
head(surv_all)
table(surv_all$fab) ## no M3

surv_aml <- surv_all %>%
  filter(fab != 'MDS normal karyotype (test set) MDS RAEB')
head(surv_aml)
surv_data <- data.frame(age = readr::parse_number(surv_aml$age),
                        time = readr::parse_number(surv_aml$time),
                        event = as.integer(str_split(surv_aml$status, ":", simplify = T)[,2]))
head(surv_data)
str(surv_data)
rownames(surv_data) <- rownames(surv_aml)
hist(as.numeric(surv_data$time))
table(surv_data$event)
## day to year
surv_data <- surv_data %>%
  filter(time >= 30)
surv_data$time <- as.numeric(surv_data$time)/365
surv_data <- surv_data %>%
  filter(time >= 30/356)
hist(surv_data$time)
head(surv_data)
str(surv_data)
## sample id
identical(colnames(dat), rownames(surv_data))
kp <- intersect(colnames(dat), rownames(surv_data))
length(kp)
surv_data <- surv_data[match(kp,rownames(surv_data)),]
head(surv_data)
dat <- dat[,match(kp,colnames(dat))]
dat[1:4, 1:4]
identical(colnames(dat), rownames(surv_data))
testset4 <- list(expression = dat,
                 survival = surv_data[,c('event','time')],
                 feature = surv_data[,c('event','time','age')])

# 2. save data----
## tidy folder
# system 里面的命令仅仅是适用于苹果电脑或者Linux电脑
# 对Windows电脑无效
system('mv *.Rdata input') 
system('mv *rda input') 

testset_list <- list(dataset1 = testset1,
                     dataset2 = testset2,
                     dataset3 = testset3,
                     dataset4 = testset4,
                     dataset5 = testset5)
save(testset, testset_list, 
     file = 'output/rdata/prepare_testset_for_model.Rdata')




