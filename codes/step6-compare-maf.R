rm(list=ls())
options(stringsAsFactors = F) 
library(maftools) 
library(dplyr)

# 0. load data----
load('output/rdata/prepare-mutation-data.Rdata')
table(dat_maf@clinical.data$Risk_level)

cg=unique(dat_maf@clinical.data$Risk_level)
cg

sub_maf_list = lapply(cg,function(x){
  kp=dat_maf@clinical.data$Risk_level == x 
  table(kp)
  cg=as.data.frame(dat_maf@clinical.data)[kp,1]
  cg
  s1.maf=read.maf(dat_maf@data[dat_maf@data$Tumor_Sample_Barcode %in% cg,],
                  clinicalData = dat_maf@clinical.data[kp,])
  s1.maf
}) 

table(dat_maf@data$Variant_Classification)
table(sub_maf_list[[1]]@data$Variant_Classification)
table(sub_maf_list[[2]]@data$Variant_Classification)
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Accent')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Silent',
  'Intron')
print(vc_cols)

names(sub_maf_list)=cg

lapply(cg,function(x){ 
  this_maf = sub_maf_list[[x]] 
  pdf(paste0('Fig4_oncoyplot_for_',x,'.pdf'))
  oncoplot(maf = this_maf , 
           colors = vc_cols,
           bgCol = "#E7EAF6FF",
           top = 20, fontSize = 0.7)
  dev.off()
}) 
library(filesstrings)

file.move('Fig4_oncoyplot_for_High.pdf', 'output/plot/')
file.move('Fig4_oncoyplot_for_Low.pdf', 'output/plot/')
