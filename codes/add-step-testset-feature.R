rm(list = ls())  
options(stringsAsFactors = F) 

library(data.table)
library(tidyverse)
library(survminer) 
library(survival)

# 0. load data----
load('output/rdata/prepare_testset_for_model.Rdata')
head(testset)
lapply(testset_list, str)

CRGs <- read.table('output/text/model_lasso_cg_genes.txt')[,1]
head(CRGs)

# 1. calculate roc----
## GSE71014 didn't upload clinical feature
add_feature <- testset_list[1:4]


# 1.1 only CRG----
CRG_list <- lapply(seq_along(add_feature), function(i){
  # i = 3
  proj <- names(add_feature)[i]
  print(proj)
  ## 1.1.1 prepare data for coxph----
  kp <- CRGs
  cdata <- add_feature[[i]]$expression[kp,]
  head(cdata)[,1:4]
  xdata <- t(cdata)
  head(xdata)[,1:4]
  ydata <- add_feature[[i]]$survival
  head(ydata)
  if(identical(rownames(xdata), rownames(ydata))){
    print('sample order is OK')}
  dat_cox <- cbind(ydata, xdata)
  head(dat_cox)
  ## prepare variale for Surv(), paste gene names
  multivariate <- paste(sort(kp), collapse = '+') 
  multivariate
  
  ## 1.1.2 multi-cox----
  s <-  paste0(' Surv(time, event) ~  ', multivariate )
  model <- coxph(as.formula(s), data = dat_cox)
  summary(model, data = dat_cox)
  risk_score <- predict(model, type = 'risk', data = dat_cox)
  dat_cox$Risk_score <- risk_score 
  
  ## 1.1.3 ROC----
  library(timeROC)
  head(dat_cox)
  new_dat <- dat_cox[, c('event', 'time','Risk_score')]
  head(new_dat)
  ## need 3 cols，time、event and risk scores
  # if (i %in% c(1,2)) {
  #   result_yrs135 <- with(new_dat, timeROC(T=time,
  #                                          delta=event,
  #                                          marker=Risk_score,
  #                                          cause=1,
  #                                          times = c(1, 3, 5),
  #                                          iid = TRUE))
  #   print(result_yrs135$AUC)
  #   return(result_yrs135)
  # }
  
#  if (i %in% c(3,4)) {
    result_yrs123 <- with(new_dat, timeROC(T=time,
                                           delta=event,
                                           marker=Risk_score,
                                           cause=1,
                                           times = c(1, 2, 3),
                                           iid = TRUE))
    print(result_yrs123$AUC)
    return(result_yrs123)
#  }
  
})
names(CRG_list) <- names(add_feature)
# 1.2 add feature----

feature_list <- lapply(seq_along(add_feature), function(i){
  # i = 3
  proj <- names(add_feature)[i]
  print(proj)
  ## 1.2.1 prepare data for coxph----
  kp <- CRGs
  cdata <- add_feature[[i]]$expression[kp,]
  head(cdata)[,1:4]
  xdata <- t(cdata)
  head(xdata)[,1:4]
  ydata <- add_feature[[i]]$feature
  head(ydata)
  if(identical(rownames(xdata), rownames(ydata))){
    print('sample order is OK')}
  dat_cox <- cbind(ydata, xdata)
  head(dat_cox)
  ## prepare variale for Surv(), paste gene names
  multivariate <- paste(c(sort(kp),'age'), collapse = '+') 
  multivariate
  
  ## 1.2.2 multi-cox----
  s <-  paste0(' Surv(time, event) ~  ', multivariate )
  model <- coxph(as.formula(s), data = dat_cox)
  summary(model, data = dat_cox)
  risk_score <- predict(model, type = 'risk', data = dat_cox)
  dat_cox$Risk_score <- risk_score 
  
  # 1.2.3 ROC----
  library(timeROC)
  head(dat_cox)
  new_dat <- dat_cox[, c('event', 'time','Risk_score')]
  head(new_dat)
  ## need 3 cols，time、event and risk scores
  # if (i %in% c(1,2)) {
  # result_yrs135 <- with(new_dat, timeROC(T=time,
  #                                 delta=event,
  #                                 marker=Risk_score,
  #                                 cause=1,
  #                                 times = c(1, 3, 5),
  #                                 iid = TRUE))
  # print(result_yrs135$AUC)
  # return(result_yrs135)
  # }

#  if (i %in% c(3,4)) {
  result_yrs123 <- with(new_dat, timeROC(T=time,
                                         delta=event,
                                         marker=Risk_score,
                                         cause=1,
                                         times = c(1, 2, 3),
                                         iid = TRUE))
  print(result_yrs123$AUC)
  return(result_yrs123)
# }
})
names(feature_list) <- names(add_feature)


# 2. table ROC----
library(reactable)
library(reactablefmtr)
head(testset)
auc_list <- lapply(seq_along(add_feature), function(n){
# n = 1
feature_auc <- data.frame(# time = paste0(c(1,3,5),' year'),
                          CRG = CRG_list[[n]]$AUC,
                          add_feature = feature_list[[n]]$AUC)
# head(feature_auc)
})
names(auc_list) <- names(add_feature)
options(digits=2)
auc_table_1 <- cbind(auc_list[[1]], auc_list[[2]]) %>% round(digits = 2)
head(auc_table_1)
colnames(auc_table_1) <- paste0(c('CRG','add_feature'),rep(c(1:2),each = 2))
head(auc_table_1)

auc_table_2 <- cbind(auc_list[[3]], auc_list[[4]]) %>% round(digits = 2)
head(auc_table_2)
colnames(auc_table_2) <- paste0(c('CRG','add_feature'),rep(c(3:4),each = 2))
head(auc_table_2)

## signifiant
dataset1_pvalue <- t.test(auc_table_1[,1],auc_table_1[,2])$`p.value`
# chisq.test(auc_table_1[,1],auc_table_1[,2])$`p.value`
dataset2_pvalue <-t.test(auc_table_1[,3],auc_table_1[,4])$`p.value`
# chisq.test(auc_table_1[,3],auc_table_1[,4])$`p.value`
dataset3_pvalue <-t.test(auc_table_2[,1],auc_table_2[,2])$`p.value`
# chisq.test(auc_table_2[,1],auc_table_2[,2])$`p.value`
dataset4_pvalue <-t.test(auc_table_2[,3],auc_table_2[,4])$`p.value`
# chisq.test(auc_table_2[,3],auc_table_2[,4])$`p.value`

auc_table <- cbind(auc_table_1,auc_table_2)
head(auc_table)
auc_table <- auc_table %>%
  add_row(add_feature1 = round(dataset1_pvalue,2),
          add_feature2 = round(dataset2_pvalue,2),
          add_feature3 = round(dataset3_pvalue,2),
          add_feature4 = round(dataset4_pvalue,2))
head(auc_table)
rownames(auc_table)[4] <- 'pvalue'

auc_reactable <- reactable(auc_table, 
          columns = list(
            CRG1 = colDef(name = "CRG"),
            add_feature1 = colDef(name = "add_feature"),
            CRG2 = colDef(name = "CRG"),
            add_feature2 = colDef(name = "add_feature"),
            CRG3 = colDef(name = "CRG"),
            add_feature3 = colDef(name = "add_feature"),
            CRG4 = colDef(name = "CRG"),
            add_feature4 = colDef(name = "add_feature")
            
          ),
          columnGroups = list(
            colGroup(name = paste0(testset$GSE[1],' | ',testset$GPL[1]), columns = c("CRG1", "add_feature1")),
            colGroup(name = paste0(testset$GSE[2],' | ',testset$GPL[2]), columns = c("CRG2", "add_feature2")),
            colGroup(name = paste0(testset$GSE[3],' | ',testset$GPL[3]), columns = c("CRG3", "add_feature3")),
            colGroup(name = paste0(testset$GSE[4],' | ',testset$GPL[4]), columns = c("CRG4", "add_feature4")))
          )

auc_reactable
# 需要调用 PhantomJS 
reactablefmtr::save_reactable_test (auc_reactable, "add_feature_auc.png")

library(filesstrings)
file.move('add_feature_auc.png', 'output/plot/', overwrite = TRUE)
library(openxlsx)
total.results <- list('Sheet1' = auc_table)
write.xlsx(total.results, file = "output/text/add_trainset_feature.xlsx")


# 3. ROC plot----
library(ggplot2)
library(colorspace)   
## prepare plot function----
hcl_palettes("qualitative", plot = TRUE)
colorspace::qualitative_hcl(5, palette = "Dynamic")
colorspace::qualitative_hcl(n = 5, palette = "Dynamic") %>% 
  colorspace::swatchplot()
plot.feature <- function(dat.compare, dat.label, tag, plot.title){
  ggplot() + 
    geom_line(data = dat.compare,aes(x = fpr, y = tpr,color = feature),size = 1) + 
    scale_color_manual(name = NULL,
                       values = c("#DB9D85", "#9DB469"),
                       labels = paste0("AUC of ",feature_factor,' is ',
                                       format(round(dat.label,2),nsmall = 2)))+
    geom_line(aes(x = c(0, 1), y = c(0,1)), 
              color = "grey",
              linetype = 'dotdash')+
    labs(x = "1 - Specificity",
         y = "Sensitivity",
         title = plot.title,
         tag = paste0(tag,' year'))+
    theme_bw()+
    theme(axis.text = element_text(size = 10, face = 'bold'),
          axis.title = element_text(size = 10, face = 'bold'),
          axis.line = element_line(linetype = 1),
          plot.title = element_text(size = 10, face = 'bold',hjust = 0.5),
          panel.grid = element_blank(),
          legend.background = element_rect(linetype = 2, size = 0.2, colour = "black"),
          legend.position = c(0.685,0.185))+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    coord_fixed()
}

feature_factor <- factor(c('CRG','add_feature'), 
                         levels = c('CRG','add_feature'))
## prepare data for plot----
data_plot <- function(CRG_list,feature_list,dn,yn){
fpr_yrs <- c(CRG_list[[dn]]$FP[,yn],feature_list[[dn]]$FP[,yn])
tpr_yrs <- c(CRG_list[[dn]]$TP[,yn],feature_list[[dn]]$TP[,yn])
auc_yrs <- c(CRG_list[[dn]]$AUC[yn],feature_list[[dn]]$AUC[yn])
dat.compare = data.frame(fpr = fpr_yrs,
                         tpr = tpr_yrs,
                         feature = rep(feature_factor, 
                                       each = nrow(CRG_list[[dn]]$TP)))
dat.label <- auc_yrs
tag <- yn
plot.title <- paste0(testset$GSE[dn],' | ',testset$GPL[dn])
plot.feature(dat.compare = dat.compare, 
             dat.label = dat.label,
             tag = tag,
             plot.title = plot.title)
}

## dataset1
plot1_1 <- data_plot(CRG_list,feature_list,dn = 1,yn =1)
plot1_1
plot1_2 <- data_plot(CRG_list,feature_list,dn = 1,yn =2)
plot1_2
plot1_3 <- data_plot(CRG_list,feature_list,dn = 1,yn =3)
plot1_3

## dataset2
plot2_1 <- data_plot(CRG_list,feature_list,dn = 2,yn =1)
plot2_1
plot2_2 <- data_plot(CRG_list,feature_list,dn = 2,yn =2)
plot2_2
plot2_3 <- data_plot(CRG_list,feature_list,dn = 2,yn =3)
plot2_3

## dataset3
plot3_1 <- data_plot(CRG_list,feature_list,dn = 3,yn =1)
plot3_1
plot3_2 <- data_plot(CRG_list,feature_list,dn = 3,yn =2)
plot3_2
plot3_3 <- data_plot(CRG_list,feature_list,dn = 3,yn =3)
plot3_3

## dataset4
plot4_1 <- data_plot(CRG_list,feature_list,dn = 4,yn =1)
plot4_1
plot4_2 <- data_plot(CRG_list,feature_list,dn = 4,yn =2)
plot4_2
plot4_3 <- data_plot(CRG_list,feature_list,dn = 4,yn =3)
plot4_3

## patch plot
library(patchwork)
add_plot_feature <- plot1_1 + plot2_1 + plot3_1 + plot4_1 + 
  plot1_2 + plot2_2 + plot3_2 + plot4_2 + 
  plot1_3 + plot2_3 + plot3_3 + plot4_3 + 
  plot_layout(ncol = 4)
add_plot_feature
dev.off()

ggsave(add_plot_feature,
       filename = file.path('output/plot/', 'add_testset_feature.pdf'),
       width = 30,
       height = 15)



