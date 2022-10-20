rm(list=ls())
options(stringsAsFactors = F) 

library(data.table)
library(dplyr)
library(survminer) 
library(survival)
# library(loose.rock)
library(futile.logger) 
library(glmSparseNet)
library(ggrisk)  
library(pheatmap) 
library(ggplot2)
library(ggsci)
library(ggstatsplot)


# 0. load data ----
load( file = 'output/rdata/prepare_data_for_lasso.Rdata' )
exprSet[1:4, 1:4]
head(phe)
identical(colnames(exprSet), phe$sample)
head(overlap_genes)

# 1. tidy data for LASSO----
## prepare xdata----
kp <- overlap_genes
cdat <- exprSet[kp, ]
xdata <- t(cdat)
xdata[1:4,1:4]

## prepare ydata----
ydata <- phe
head(ydata)
identical(rownames(xdata), ydata$sample)

# 2. run lasso model----
source('codes/run_LASSO_cox.R') 

lasso_cox(xdata,
          ydata,savDir = 'output')









