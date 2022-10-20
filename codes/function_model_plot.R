
## 2.1 forest----
pancancer_forest <- function(model, dat_cox, savDir='./'){
library(survminer)
options(scipen=1)
ggforest(model, data = dat_cox, 
         main = "Hazard ratio", 
         cpositions = c(0.06, 0.22, 0.4), 
         fontsize = 1.0, 
         refLabel = "1", noDigits = 4)
ggsave(file.path(savDir, 'multicox_forest.pdf'), height = 5, width = 13)
}

## dotchart----
pancancer_dotchart <- function(df,savDir='./'){
  options(scipen=1)
  df$logHR=log(df$HR)
  df$Sig = ifelse(df$p.val < 0.05,'*','no')
  table(df$Sig)
  df$cancer=rownames(df)
  ggdotchart(df, y = "logHR", x = "cancer",
             color = "Sig",       
             sorting = 'descending',
             group = 'p.val',
             add = "segments", 
             palette = "nejm",
             rotate = TRUE,
             dot.size = 4,
             xlab = "",
             ylab = 'logHR',
             ggtheme = theme_bw())+ 
    geom_vline(xintercept= 0,lty=4,col="grey",lwd=0.8) +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 12, face = 'bold'),
          axis.text = element_text(size = 10, face = 'bold'),
          legend.position = 'top')
ggsave(file.path(savDir, 'multicox_dotchart.pdf'), height = 4, width = 6)
}

## 2.2 nomo----
pancancer_nomo <- function(dat_cox,savDir='./'){
  attach(dat_cox)
library(rms)
dc <- datadist(dat_cox);dc
options(datadist="dc")
boxplot(dat_cox$time)
# Cox Proportional Hazards Model and Extensions , cph {rms}
f <- cph(as.formula(s),
         x=T, y=T, surv=T, 
         data=dat_cox, time.inc = 15)
summary(f) 

# Fit Proportional Hazards Regression Model , coxph {survival} 
surv<- Survival(f)
nom <- nomogram(f, fun=list(function(x) surv(1, x), 
                            function(x) surv(3, x),
                            function(x) surv(5, x)), 
                lp=F, funlabel=c("1-year survival", 
                                 "3-yearsurvival",
                                 "5-year survival"))
# maxscale=10, 
# fun.at=c(0.95,0.9, 0.85, 0.8, 0.75, 0.7, 0.6, 0.5)
pdf(file.path(savDir, 'multicox_nomogram.pdf'), width = 15)
plot(nom)
dev.off()
}

## 2.3 ROC----
pancancer_roc_yrs135 <- function(new_dat,savDir='./'){
library(timeROC)
library(ggstatsplot)
## need 3 cols，time、event and risk scores
result <- with(new_dat, timeROC(T=time,
                                delta=event,
                                marker=Risk_score,
                                cause=1,
                                times = c(1, 3, 5),
                                iid = TRUE))
#identical(c(result$TP[,1],result$TP[,2],result$TP[,3]),as.numeric(result$TP))
dat = data.frame(fpr = as.numeric(result$FP),
                 tpr = as.numeric(result$TP),
                 time = rep(as.factor(c(1,3,5)),each = nrow(result$TP)))

## 将geom_line()改为geom_smooth(method = "loess")
## 在数学中有种算法叫“样条插补法”，这种方法可以获得过点的平滑曲线
ggplot() + 
  geom_smooth(data = dat, 
              aes(x = fpr, y = tpr, color = time), 
              size = 1,
              method = "loess",
              se = FALSE) + 
  scale_color_manual(name = NULL,
                     values = c("#92C5DE", "#F4A582", "#66C2A5"),
                     labels = paste0("AUC of ",c(1,3,5),"-year survival: ",
                                     format(round(result$AUC,2),nsmall = 2)))+
  geom_line(aes(x = c(0, 1), y = c(0,1)), 
            color = "grey",
            linetype = 'dotdash')+
  theme_bw()+
  theme(axis.text = element_text(size = 10, face = 'bold'),
        axis.title = element_text(size = 10, face = 'bold'),
        axis.line = element_line(linetype = 1),
        panel.grid = element_blank(),
        legend.background = element_rect(linetype = 2, size = 0.2, colour = "black"),
        legend.position = c(0.665,0.135))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  coord_fixed()
ggsave(file.path(savDir, 'smooth_ROC.pdf'))
}
pancancer_roc_yrs123 <- function(new_dat,savDir='./'){
  library(timeROC)
  library(ggstatsplot)
  ## need 3 cols，time、event and risk scores
  result <- with(new_dat, timeROC(T=time,
                                  delta=event,
                                  marker=Risk_score,
                                  cause=1,
                                  times = c(1, 2, 3),
                                  iid = TRUE))
  #identical(c(result$TP[,1],result$TP[,2],result$TP[,3]),as.numeric(result$TP))
  dat = data.frame(fpr = as.numeric(result$FP),
                   tpr = as.numeric(result$TP),
                   time = rep(as.factor(c(1,3,5)),each = nrow(result$TP)))
  
  ## 将geom_line()改为geom_smooth(method = "loess")
  ## 在数学中有种算法叫“样条插补法”，这种方法可以获得过点的平滑曲线
  ggplot() + 
    geom_smooth(data = dat, 
                aes(x = fpr, y = tpr, color = time), 
                size = 1,
                method = "loess",
                se = FALSE) + 
    scale_color_manual(name = NULL,
                       values = c("#92C5DE", "#F4A582", "#66C2A5"),
                       labels = paste0("AUC of ",c(1,2,3),"-year survival: ",
                                       format(round(result$AUC,2),nsmall = 2)))+
    geom_line(aes(x = c(0, 1), y = c(0,1)), 
              color = "grey",
              linetype = 'dotdash')+
    theme_bw()+
    theme(axis.text = element_text(size = 10, face = 'bold'),
          axis.title = element_text(size = 10, face = 'bold'),
          axis.line = element_line(linetype = 1),
          panel.grid = element_blank(),
          legend.background = element_rect(linetype = 2, size = 0.2, colour = "black"),
          legend.position = c(0.665,0.135))+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    labs(x = "1 - Specificity",
         y = "Sensitivity")+
    coord_fixed()
  ggsave(file.path(savDir, 'smooth_ROC.pdf'))
}
## 2.4 km----
pancancer_km <- function(new_dat, savDir='./'){
  attach(new_dat)

sfit <- survfit(Surv(time, event)~Risk_level, data=new_dat)
sfit
summary(sfit)

## more complicate figures.
survp=ggsurvplot(
  sfit,                     # survfit object with calculated statistics.
  legend.title = 'Risk level', 
  # legend = "top",#图例位置
  # legend.labs = c('High', 'Low'),
  pval = T, #在图上添加log rank检验的p值
  risk.table = TRUE, 
  risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
  xlab = "Time in years", #x轴标题
  # xlim = c(0, 10), #展示x轴的范围
  break.time.by = 1, #x轴间隔
  size = 1.5, #线条大小
  ggtheme = theme_bw(),
  palette="nejm", #配色
)
print(survp)
pdf(file.path(savDir, 'multicox_KM.pdf'), onefile = F, width = 10)
print(survp)
dev.off()
}

## 2.5 ggrisk----
pancancer_ggrisk <- function(model,savDir='./'){
## https://cran.r-project.org/web/packages/ggrisk/ggrisk.pdf
## https://cloud.tencent.com/developer/article/1765625
library(ggrisk)
# save in pdf
pdf(file.path(savDir, 'riskscore.pdf'), onefile = F)
ggrisk(model,
       color.A = c(low = "#0B5E9D", high = "#EA5205"),
       color.B = c(code.0 = "#0B5E9D", code.1 = "#EA5205"),
       color.C = c(low = "#0B5E9D", median = "white", high = "#EA5205"))
dev.off() 
}



