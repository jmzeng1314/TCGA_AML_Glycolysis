rm(list=ls())
options(stringsAsFactors = F)

library(survminer) 
library(survival)
library(timeROC)
library(rms)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggstatsplot)

# 0. load data----
load('output/rdata/coefs.v_lasso_model.Rdata')
xdata[1:4, 1:4]
head(ydata)
head(names(coefs.v))

# 1. tidy data for multicox----
if(identical(rownames(xdata), ydata$sample)){
  print('sample order is OK')}
kp <- names(coefs.v)
dat_cox <- cbind(ydata[,-4], xdata[,kp])
head(dat_cox)

if(F){
  ## z-score
  cg <- sort(names(coefs.v));cg
  n <- apply(xdata,2,scale)[,cg] %>%
    as.data.frame()
  head(n)
  rownames(n) <- rownames(xdata)
  boxplot(n)
  #colnames(n)= gsub('-','_',colnames(n))
  dat_cox <- cbind(ydata[,-4],n)
  head(dat_cox)
}

## prepare variale for Surv(), paste gene names
multivariate <- paste(sort(kp), collapse = '+') 
multivariate

## 'model' is for ggforest and ggrisk----
s <-  paste0(' Surv(time, event) ~  ', multivariate )
model <- coxph(as.formula(s), data = dat_cox)
# ggrisk::ggrisk(model)

## 'df' is for dotchart----
HR_result <-  summary(model, data = dat_cox)
df <- HR_result$coefficients %>%
  as.data.frame() %>%
  select(2, 5)
head(df)
colnames(df) <- c('HR', 'p.val')

## 'dc' is for ggrisk and nomo---- 
dc <- datadist(dat_cox)

## 'risk_score' is for roc----
# predict dataself
risk_score <- predict(model, type = 'risk', data = dat_cox)
new_dat <- dat_cox[, c('event', 'time')]
new_dat$Risk_score <- risk_score 
head(new_dat)

## 'Risk_level' is for  km----
new_dat$Risk_level <- as.factor(ifelse(new_dat$Risk_score > median(new_dat$Risk_score),'High','Low'))
head(new_dat)

summary(model)

dir.create('output/plot/trainset',recursive = T)
source('codes/function_model_plot.R')
proj <- 'trainset'
pancancer_forest(model, dat_cox, savDir= paste0('output/plot/',proj))
pancancer_nomo(dat_cox,savDir= paste0('output/plot/',proj))
pancancer_roc(new_dat, savDir= paste0('output/plot/',proj))
pancancer_km(new_dat,savDir= paste0('output/plot/',proj))
pancancer_dotchart(df, savDir= paste0('output/plot/',proj))

# pancancer_ggrisk(model, savDir= paste0('output/plot/',proj))



