rm(list=ls())
options(stringsAsFactors = F)

library(tidyverse)
library(stringr)
library(data.table)
library(ggplot2)
library(ggsci)
library(ggstatsplot)

# 0. load data----
load('output/rdata/prepare_testset_for_model.Rdata')
head(testset)
lapply(testset_list, str)
signature <- read.table('output/text/model_lasso_cg_genes.txt')[,1]
head(signature)
dir.create('output/plot/testset',recursive = T)

# 1. do Model----
library(survminer) 
library(survival)
library(timeROC)
library(rms)

lapply(seq_along(testset_list), function(i){
# lapply(c(4:5), function(i){
# i = 3
proj <- names(testset_list)[i]
print(proj)
dir.create(paste0('output/plot/testset/',proj),recursive = T)

# 1.1 prepare data for coxph----
kp <- signature
cdata <- testset_list[[i]]$expression[kp,]
head(cdata)[,1:4]
xdata <- t(cdata)
head(xdata)[,1:4]
ydata <- testset_list[[i]]$survival
head(ydata)
if(identical(rownames(xdata), rownames(ydata))){
  print('sample order is OK')}
dat_cox <- cbind(ydata, xdata)
head(dat_cox)

## prepare variale for Surv(), paste gene names
multivariate <- paste(sort(kp), collapse = '+') 
multivariate

## 1.2 'model' is for ggforest and ggrisk----
s <-  paste0(' Surv(time, event) ~  ', multivariate )
model <- coxph(as.formula(s), data = dat_cox)

## 1.3 'df' is for dotchart----
HR_result <-  summary(model, data = dat_cox)
df <- HR_result$coefficients %>%
  as.data.frame() %>%
  select(2, 5)
head(df)
colnames(df) <- c('HR', 'p.val')

## 1.4 'dc' is for ggrisk and nomo---- 
dc <- datadist(dat_cox)

## 1.5 'risk_score' is for roc----
# predict dataself
risk_score <- predict(model, type = 'risk', data = dat_cox)
new_dat <- dat_cox[, c('event', 'time')]
new_dat$Risk_score <- risk_score 
head(new_dat)

## 1.6 'Risk_level' is for  km----
new_dat$Risk_level <- as.factor(ifelse(new_dat$Risk_score > median(new_dat$Risk_score),'High','Low'))
head(new_dat)

# 2. plot----
source('codes/function_model_plot.R')
pancancer_forest(model, dat_cox, savDir= paste0('output/plot/testset/',proj))
# pancancer_nomo(dat_cox,savDir= paste0('output/plot/testset/',proj))

if (i %in% c(3,4) ) {
  pancancer_roc_yrs123(new_dat, savDir= paste0('output/plot/testset/',proj))
}

if (i %in% c(1,2,5) ) {
  pancancer_roc_yrs135(new_dat, savDir= paste0('output/plot/testset/',proj))
}

pancancer_km(new_dat,savDir= paste0('output/plot/testset/',proj))
pancancer_dotchart(df, savDir= paste0('output/plot/testset/',proj))
# pancancer_ggrisk(model, savDir= paste0('output/plot/testset/',proj))
# return(new_dat)
})





