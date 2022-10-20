rm(list = ls())  
options(stringsAsFactors = F) 

library(data.table)
library(tidyverse)
library(survminer) 
library(survival)

# 0. load data----
load( file = 'output/rdata/prepare_data_for_lasso.Rdata' )
exprSet[1:4, 1:4]
head(phe)
identical(colnames(exprSet), phe$sample)

CRGs <- read.table('output/text/model_lasso_cg_genes.txt')[,1]
head(CRGs)


# 1. check other signatures in expression----
Chen_2020 <- c('CASP3', 'CHAF1B', 'KLHL24', 'OPTN', 'VEGFA', 'VPS37C')
Qu_2020 <- c('FLT3', 'CD177', 'TTPAL')
Cai_2020 <- c('EEF1A1','RPLP2','RPL19')
Jiang_2021 <- c('STAT1', 'BATF', 'EML4')
# Liu_2022 <- c('ELANE','CASP3','CASP6','GPX4','CASP1','CASP9','AIM2','PYCARD')
# Li_2022 <- c('LDLRAP1','PNPLA6','DGKA','PLA2G4A','CBR1','EBP')

length(intersect(Chen_2020, rownames(exprSet)))
length(intersect(Qu_2020, rownames(exprSet)))
length(intersect(Cai_2020, rownames(exprSet)))
length(intersect(Jiang_2021, rownames(exprSet)))
#length(intersect(Liu_2022, rownames(exprSet)))
# length(intersect(Li_2022, rownames(exprSet)))
comparative <- list(CRGs = CRGs,
                    Chen_2020 = Chen_2020,
                    Qu_2020 = Qu_2020,
                    Cai_2020 = Cai_2020,
                    Jiang_2021 = Jiang_2021)
# 1. multivarible-cox----
## 1.1 tidy data for multicox----
## prepare xdata----
comparative_list <- lapply(seq_along(comparative), function(i){

# i = 1
  print(names(comparative)[i])
kp <- comparative[[i]]
cdat <- exprSet[kp, ]
xdata <- t(cdat)
head(xdata)

## prepare ydata----
ydata <- phe
head(ydata)
identical(rownames(xdata), ydata$sample)

dat_cox <- cbind(ydata,xdata)
head(dat_cox)

## prepare variale for Surv(), paste gene names
multivariate <- paste(sort(colnames(xdata)), collapse = '+') 

## 1.2 multi-cox----
attach(dat_cox)
s <-  paste0(' Surv(time, event) ~  ', multivariate )
model <- coxph(as.formula(s), data = dat_cox )
summary(model, data = dat_cox)
risk_score <- predict(model, type = 'risk', data = dat_cox)
dat_cox$Risk_score <- risk_score 

# 2. ROC----
library(timeROC)
## 2.1 calculate ROC----
head(dat_cox)
new_dat <- dat_cox[, c('event', 'time','Risk_score')]
head(new_dat)
## need 3 cols，time、event and risk scores
result <- with(new_dat, timeROC(T=time,
                                delta=event,
                                marker=Risk_score,
                                cause=1,
                                times = c(1, 3, 5),
                                iid = TRUE))
print(result$AUC)
return(result)
})
names(comparative_list) <- names(comparative)


## 2.2 comparative ROC plot----
library(ggplot2)
library(colorspace)   
## prepare plot function----
hcl_palettes("qualitative", plot = TRUE)
colorspace::qualitative_hcl(5, palette = "Dynamic")
colorspace::qualitative_hcl(n = 5, palette = "Dynamic") %>% 
  colorspace::swatchplot()
plot.comparative <- function(dat.compare, dat.label, plot.title){
  ggplot() + 
    geom_line(data = dat.compare,aes(x = fpr, y = tpr,color = signature),size = 1) + 
    scale_color_manual(name = NULL,
                       values = c("red", "#9DB469", "#3DBEAB",
                                  '#87AEDF','#DA95CC'),
                       labels = paste0("AUC of ",names(comparative),' is ',
                                       format(round(dat.label,2),nsmall = 2)))+
    geom_line(aes(x = c(0, 1), y = c(0,1)), 
              color = "grey",
              linetype = 'dotdash')+
    labs(x = "1 - Specificity",
         y = "Sensitivity",
         title = plot.title)+
    theme_bw()+
    theme(axis.text = element_text(size = 10, face = 'bold'),
          axis.title = element_text(size = 10, face = 'bold'),
          axis.line = element_line(linetype = 1),
          plot.title = element_text(size = 10, face = 'bold',hjust = 0.5),
          panel.grid = element_blank(),
          legend.background = element_rect(linetype = 2, size = 0.2, colour = "black"),
          legend.position = c(0.785,0.185))+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    coord_fixed()
}


## prepare data for plot----
fpr_yrs1 <- unlist(lapply(seq_along(comparative_list), function(a){
  # a = 1
  comparative_list[[a]]$FP[,1]
}))

fpr_yrs3 <- unlist(lapply(seq_along(comparative_list), function(a){
  # a = 1
  comparative_list[[a]]$FP[,2]
}))
fpr_yrs5 <- unlist(lapply(seq_along(comparative_list), function(a){
  # a = 1
  comparative_list[[a]]$FP[,3]
}))

tpr_yrs1 <- unlist(lapply(seq_along(comparative_list), function(b){
  # b = 1
  comparative_list[[b]]$TP[,1]
}))
tpr_yrs3 <- unlist(lapply(seq_along(comparative_list), function(b){
  # b = 1
  comparative_list[[b]]$TP[,2]
}))
tpr_yrs5 <- unlist(lapply(seq_along(comparative_list), function(b){
  # b = 1
  comparative_list[[b]]$TP[,3]
}))

auc_yrs1 <- unlist(lapply(seq_along(comparative_list), function(d){
  # d = 1
  comparative_list[[d]]$AUC[1]
}))
auc_yrs3 <- unlist(lapply(seq_along(comparative_list), function(d){
  # d = 1
  comparative_list[[d]]$AUC[2]
}))
auc_yrs5 <- unlist(lapply(seq_along(comparative_list), function(d){
  # d = 1
  comparative_list[[d]]$AUC[3]
}))

signature_factor <- factor(names(comparative), 
                           levels = names(comparative))


## 1 year
dat.compare = data.frame(fpr = fpr_yrs1,
                         tpr = tpr_yrs1,
                         signature = rep(signature_factor, 
                                         each = nrow(comparative_list[[1]]$TP))
                         )
dat.label <- auc_yrs1
plot.title <- 'ROC at 1 year'
plot.yrs1 <- plot.comparative(dat.compare = dat.compare, 
                              dat.label = dat.label,
                              plot.title = plot.title)
plot.yrs1
## 3 year
dat.compare = data.frame(fpr = fpr_yrs3,
                         tpr = tpr_yrs3,
                         signature = rep(signature_factor, 
                                         each = nrow(comparative_list[[1]]$TP))
)
plot.title <- 'ROC at 3 year'
dat.label <- auc_yrs3
plot.yrs3 <- plot.comparative(dat.compare = dat.compare, 
                              dat.label = dat.label,
                              plot.title = plot.title)
plot.yrs3

## 5 year
dat.compare = data.frame(fpr = fpr_yrs5,
                         tpr = tpr_yrs5,
                         signature = rep(signature_factor, 
                                         each = nrow(comparative_list[[1]]$TP))
)
plot.title <- 'ROC at 5 year'
dat.label <- auc_yrs5
plot.yrs5 <- plot.comparative(dat.compare = dat.compare, 
                              dat.label = dat.label,
                              plot.title = plot.title)
plot.yrs5

## patch plot
library(patchwork)
add_plot_comparative <- plot.yrs1+plot.yrs3+plot.yrs5
add_plot_comparative
ggsave(add_plot_comparative,
       filename = file.path('output/plot/', 'add_signature_comparative.pdf'),
       width = 18)



