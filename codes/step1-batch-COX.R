rm(list = ls())  
options(stringsAsFactors = F) 

library(survminer) 
library(survival)
library(tidyverse)

# 0. load data ----
load(file = 'output/rdata/prepare_data_for_cox.Rdata')
mat[1:4, 1:4]
head(surv_data)
identical(colnames(mat), surv_data$sample)

# 1. prepare data for coxph----
exprSet = log2(edgeR::cpm(mat)+1)
dim(exprSet)
exprSet[1:4,1:4]

phe <- surv_data
head(phe) 

## 1.1 use coxph----
## http://www.sthda.com/english/wiki/cox-proportional-hazards-model
mySurv <- with(phe, Surv(time, event))

cox_results <-apply(exprSet , 1 , function(gene){
  # gene= as.numeric(exprSet[1,])
  group=ifelse(gene>median(gene),'high','low') 
  if( length(table(group))<2)
    return(NULL)
  survival_dat <- data.frame(group=group,# stage=phe$stage,
                             stringsAsFactors = F)
  m=coxph(mySurv ~ group, 
          data =  survival_dat)
  
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
  return(tmp['grouplow',])
  
})

# 2. specify the value----
cox_results=t(cox_results)
head(cox_results)
table(cox_results[,4]<0.01)
table(cox_results[,4]<0.05)

## overlap with target genes
p_threshold <- 0.01
cox_threshold <- cox_results %>% as.data.frame() %>%
  filter(p < p_threshold)
head(cox_threshold)
dim(cox_threshold)
overlap_genes <- intersect(genes_target, rownames(cox_threshold))
length(overlap_genes)
overlap_genes

# 3.save data ----
save(overlap_genes, exprSet, phe,
     file = 'output/rdata/prepare_data_for_lasso.Rdata')










 