
#this file is created to run more real datasets

library(survAUC)
library(reshape2)
library(ggplot2)
library(readxl)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(mstate)
library(MASS)
library(survival)
library(MatchIt)
library(mice)
library(DT)
library(EnvStats)
library(pbmcapply)
library(FNN)

library(pec)
library(MultiAssayExperiment)
library(TCGAbiolinks)

library(survivalROC)
library(survival)
#library(CoxBoost) #archrived package
library(limma)
library(survivalsvm)
library(survminer)
library(randomForestSRC)
library(glmnet)
library(rms)
library(penalized)
library(riskRegression)
library(pROC)
library(ROCR)
library(cvTools)
library(parallel)
library(MTLR)
library(icenReg)#was trying to get data (miceData), however, there is no covariates with that interval censored data
library(mcsurvdata)
library(ExperimentHub)
library(hdnom)
library(Hmisc)
library(ggrepel)

source("functions.R")

# this is a server run script which takes the codes from code13 to run only
stability_selection_fun1=function(i,data,ex_validation,timess){
  set.seed(i)
  updated_data=bootstrap_fun(data,i)
  censoring_rate_train=table(updated_data[[1]]$status)[1]/(table(updated_data[[1]]$status)[1]+table(updated_data[[1]]$status)[2])
  censoring_rate_test=table(updated_data[[2]]$status)[1]/(table(updated_data[[2]]$status)[1]+table(updated_data[[2]]$status)[2])
  rr=lasso_feature_selection_withoutcal(updated_data,ex_validation,timess)
  coeff=rr[[1]]
  cindex=rr[[2]]
  length_coef=rr[[3]]
  bs=rr[[4]]
  auc=rr[[5]]
  distance_train=ks.test(updated_data[[1]]$time,data$time)
  distance_test=ks.test(updated_data[[2]]$time,data$time)
  lambda=rr[[16]]
  npasses=rr[[17]]
  devratio=rr[[18]]
  cindex_train=rr[[19]]
  bs_train=rr[[20]]
  cindex_ex=rr[[21]]
  bs_ex=rr[[22]]
  auc_ex=rr[[23]]
  auc_ex_all=c(unlist(rr[[24]]),unlist(rr[[25]]),unlist(rr[[26]]),unlist(rr[[27]]),unlist(rr[[28]]),unlist(rr[[29]]),unlist(rr[[30]]),unlist(rr[[31]]),unlist(rr[[32]]),unlist(rr[[33]]))
  auc_te_all=c(unlist(rr[[6]]),unlist(rr[[7]]),unlist(rr[[8]]),unlist(rr[[9]]),unlist(rr[[10]]),unlist(rr[[11]]),unlist(rr[[12]]),unlist(rr[[13]]),unlist(rr[[14]]),unlist(rr[[15]]))
  # cal_train=rr[[34]]
  # cal_test=rr[[35]]
  # cal_ex=rr[[36]]
  # return(list(censoring_rate_train,coeff,cindex,length_coef,bs,auc,censoring_rate_test,distance_train,distance_test,lambda,npasses,devratio,cindex_train,bs_train,cindex_ex,bs_ex,auc_ex,auc_ex_all,auc_te_all,cal_train,cal_test,cal_ex))
  return(list(censoring_rate_train,coeff,cindex,length_coef,bs,auc,censoring_rate_test,distance_train,distance_test,lambda,npasses,devratio,cindex_train,bs_train,cindex_ex,bs_ex,auc_ex,auc_ex_all,auc_te_all))
}

#write into one function to run in parallell
stability_selection_fun2=function(i,data,ex_validation,timess,prob){
  set.seed(i)
  updated_data=bootstrap_fun2(data,i,desired_proportion_A=prob)
  censoring_rate_train=table(updated_data[[1]]$status)[1]/(table(updated_data[[1]]$status)[1]+table(updated_data[[1]]$status)[2])
  censoring_rate_test=table(updated_data[[2]]$status)[1]/(table(updated_data[[2]]$status)[1]+table(updated_data[[2]]$status)[2])
  rr=lasso_feature_selection_withoutcal(updated_data,ex_validation,timess)
  coeff=rr[[1]]
  cindex=rr[[2]]
  length_coef=rr[[3]]
  bs=rr[[4]]
  auc=rr[[5]]
  distance_train=ks.test(updated_data[[1]]$time,data$time)
  distance_test=ks.test(updated_data[[2]]$time,data$time)
  lambda=rr[[16]]
  npasses=rr[[17]]
  devratio=rr[[18]]
  cindex_train=rr[[19]]
  bs_train=rr[[20]]
  cindex_ex=rr[[21]]
  bs_ex=rr[[22]]
  auc_ex=rr[[23]]
  auc_ex_all=c(unlist(rr[[24]]),unlist(rr[[25]]),unlist(rr[[26]]),unlist(rr[[27]]),unlist(rr[[28]]),unlist(rr[[29]]),unlist(rr[[30]]),unlist(rr[[31]]),unlist(rr[[32]]),unlist(rr[[33]]))
  auc_te_all=c(unlist(rr[[6]]),unlist(rr[[7]]),unlist(rr[[8]]),unlist(rr[[9]]),unlist(rr[[10]]),unlist(rr[[11]]),unlist(rr[[12]]),unlist(rr[[13]]),unlist(rr[[14]]),unlist(rr[[15]]))
  # cal_train=rr[[34]]
  # cal_test=rr[[35]]
  # cal_ex=rr[[36]]
  # return(list(censoring_rate_train,coeff,cindex,length_coef,bs,auc,censoring_rate_test,distance_train,distance_test,lambda,npasses,devratio,cindex_train,bs_train,cindex_ex,bs_ex,auc_ex,auc_ex_all,auc_te_all,cal_train,cal_test,cal_ex))
  return(list(censoring_rate_train,coeff,cindex,length_coef,bs,auc,censoring_rate_test,distance_train,distance_test,lambda,npasses,devratio,cindex_train,bs_train,cindex_ex,bs_ex,auc_ex,auc_ex_all,auc_te_all))
}

#write into one function to run in parallell
stability_selection_fun3=purrr::possibly(function(i,data,ex_validation,timess,m_percentage){
  set.seed(i)
  updated_data=moutofn_bootstrap_fun4(data,i,m_percentage)
  censoring_rate_train=table(updated_data[[1]]$status)[1]/(table(updated_data[[1]]$status)[1]+table(updated_data[[1]]$status)[2])
  censoring_rate_test=table(updated_data[[2]]$status)[1]/(table(updated_data[[2]]$status)[1]+table(updated_data[[2]]$status)[2])
  rr=lasso_feature_selection_withoutcal(updated_data,ex_validation,timess)
  coeff=rr[[1]]
  cindex=rr[[2]]
  length_coef=rr[[3]]
  bs=rr[[4]]
  auc=rr[[5]]
  distance_train=ks.test(updated_data[[1]]$time,data$time)
  distance_test=ks.test(updated_data[[2]]$time,data$time)
  lambda=rr[[16]]
  npasses=rr[[17]]
  devratio=rr[[18]]
  cindex_train=rr[[19]]
  bs_train=rr[[20]]
  cindex_ex=rr[[21]]
  bs_ex=rr[[22]]
  auc_ex=rr[[23]]
  auc_ex_all=c(unlist(rr[[24]]),unlist(rr[[25]]),unlist(rr[[26]]),unlist(rr[[27]]),unlist(rr[[28]]),unlist(rr[[29]]),unlist(rr[[30]]),unlist(rr[[31]]),unlist(rr[[32]]),unlist(rr[[33]]))
  auc_te_all=c(unlist(rr[[6]]),unlist(rr[[7]]),unlist(rr[[8]]),unlist(rr[[9]]),unlist(rr[[10]]),unlist(rr[[11]]),unlist(rr[[12]]),unlist(rr[[13]]),unlist(rr[[14]]),unlist(rr[[15]]))
  # cal_train=rr[[34]]
  #  cal_test=rr[[35]]
  #  cal_ex=rr[[36]]
  #  return(list(censoring_rate_train,coeff,cindex,length_coef,bs,auc,censoring_rate_test,distance_train,distance_test,lambda,npasses,devratio,cindex_train,bs_train,cindex_ex,bs_ex,auc_ex,auc_ex_all,auc_te_all,cal_train,cal_test,cal_ex))
  return(list(censoring_rate_train,coeff,cindex,length_coef,bs,auc,censoring_rate_test,distance_train,distance_test,lambda,npasses,devratio,cindex_train,bs_train,cindex_ex,bs_ex,auc_ex,auc_ex_all,auc_te_all))
},otherwise = NA)

stability_selection_fun4=purrr::possibly(function(i,data,ex_validation,timess,m_percentage){
  set.seed(i)
  updated_data=moutofn_bootstrap_fun_wr4(data,i,m_percentage)
  censoring_rate_train=table(updated_data[[1]]$status)[1]/(table(updated_data[[1]]$status)[1]+table(updated_data[[1]]$status)[2])
  censoring_rate_test=table(updated_data[[2]]$status)[1]/(table(updated_data[[2]]$status)[1]+table(updated_data[[2]]$status)[2])
  rr=lasso_feature_selection_withoutcal(updated_data,ex_validation,timess)
  coeff=rr[[1]]
  cindex=rr[[2]]
  length_coef=rr[[3]]
  bs=rr[[4]]
  auc=rr[[5]]
  distance_train=ks.test(updated_data[[1]]$time,data$time)
  distance_test=ks.test(updated_data[[2]]$time,data$time)
  lambda=rr[[16]]
  npasses=rr[[17]]
  devratio=rr[[18]]
  cindex_train=rr[[19]]
  bs_train=rr[[20]]
  cindex_ex=rr[[21]]
  bs_ex=rr[[22]]
  auc_ex=rr[[23]]
  auc_ex_all=c(unlist(rr[[24]]),unlist(rr[[25]]),unlist(rr[[26]]),unlist(rr[[27]]),unlist(rr[[28]]),unlist(rr[[29]]),unlist(rr[[30]]),unlist(rr[[31]]),unlist(rr[[32]]),unlist(rr[[33]]))
  auc_te_all=c(unlist(rr[[6]]),unlist(rr[[7]]),unlist(rr[[8]]),unlist(rr[[9]]),unlist(rr[[10]]),unlist(rr[[11]]),unlist(rr[[12]]),unlist(rr[[13]]),unlist(rr[[14]]),unlist(rr[[15]]))
  # cal_train=rr[[34]]
  #  cal_test=rr[[35]]
  #  cal_ex=rr[[36]]
  #  return(list(censoring_rate_train,coeff,cindex,length_coef,bs,auc,censoring_rate_test,distance_train,distance_test,lambda,npasses,devratio,cindex_train,bs_train,cindex_ex,bs_ex,auc_ex,auc_ex_all,auc_te_all,cal_train,cal_test,cal_ex))
  return(list(censoring_rate_train,coeff,cindex,length_coef,bs,auc,censoring_rate_test,distance_train,distance_test,lambda,npasses,devratio,cindex_train,bs_train,cindex_ex,bs_ex,auc_ex,auc_ex_all,auc_te_all))
},otherwise = NA)

#without restriction
#write into one function to run in parallell
stability_selection_fun5=purrr::possibly(function(i,data,ex_validation,timess,m_percentage){
  set.seed(i)
  updated_data=moutofn_bootstrap_fun(data,i,m_percentage)
  censoring_rate_train=table(updated_data[[1]]$status)[1]/(table(updated_data[[1]]$status)[1]+table(updated_data[[1]]$status)[2])
  censoring_rate_test=table(updated_data[[2]]$status)[1]/(table(updated_data[[2]]$status)[1]+table(updated_data[[2]]$status)[2])
  rr=lasso_feature_selection_withoutcal(updated_data,ex_validation,timess)
  coeff=rr[[1]]
  cindex=rr[[2]]
  length_coef=rr[[3]]
  bs=rr[[4]]
  auc=rr[[5]]
  distance_train=ks.test(updated_data[[1]]$time,data$time)
  distance_test=ks.test(updated_data[[2]]$time,data$time)
  lambda=rr[[16]]
  npasses=rr[[17]]
  devratio=rr[[18]]
  cindex_train=rr[[19]]
  bs_train=rr[[20]]
  cindex_ex=rr[[21]]
  bs_ex=rr[[22]]
  auc_ex=rr[[23]]
  auc_ex_all=c(unlist(rr[[24]]),unlist(rr[[25]]),unlist(rr[[26]]),unlist(rr[[27]]),unlist(rr[[28]]),unlist(rr[[29]]),unlist(rr[[30]]),unlist(rr[[31]]),unlist(rr[[32]]),unlist(rr[[33]]))
  auc_te_all=c(unlist(rr[[6]]),unlist(rr[[7]]),unlist(rr[[8]]),unlist(rr[[9]]),unlist(rr[[10]]),unlist(rr[[11]]),unlist(rr[[12]]),unlist(rr[[13]]),unlist(rr[[14]]),unlist(rr[[15]]))
  # cal_train=rr[[34]]
  #  cal_test=rr[[35]]
  #  cal_ex=rr[[36]]
  #  return(list(censoring_rate_train,coeff,cindex,length_coef,bs,auc,censoring_rate_test,distance_train,distance_test,lambda,npasses,devratio,cindex_train,bs_train,cindex_ex,bs_ex,auc_ex,auc_ex_all,auc_te_all,cal_train,cal_test,cal_ex))
  return(list(censoring_rate_train,coeff,cindex,length_coef,bs,auc,censoring_rate_test,distance_train,distance_test,lambda,npasses,devratio,cindex_train,bs_train,cindex_ex,bs_ex,auc_ex,auc_ex_all,auc_te_all))
},otherwise = NA)

stability_selection_fun6=purrr::possibly(function(i,data,ex_validation,timess,m_percentage){
  set.seed(i)
  updated_data=moutofn_bootstrap_fun_wr(data,i,m_percentage)
  censoring_rate_train=table(updated_data[[1]]$status)[1]/(table(updated_data[[1]]$status)[1]+table(updated_data[[1]]$status)[2])
  censoring_rate_test=table(updated_data[[2]]$status)[1]/(table(updated_data[[2]]$status)[1]+table(updated_data[[2]]$status)[2])
  rr=lasso_feature_selection_withoutcal(updated_data,ex_validation,timess)
  coeff=rr[[1]]
  cindex=rr[[2]]
  length_coef=rr[[3]]
  bs=rr[[4]]
  auc=rr[[5]]
  distance_train=ks.test(updated_data[[1]]$time,data$time)
  distance_test=ks.test(updated_data[[2]]$time,data$time)
  lambda=rr[[16]]
  npasses=rr[[17]]
  devratio=rr[[18]]
  cindex_train=rr[[19]]
  bs_train=rr[[20]]
  cindex_ex=rr[[21]]
  bs_ex=rr[[22]]
  auc_ex=rr[[23]]
  auc_ex_all=c(unlist(rr[[24]]),unlist(rr[[25]]),unlist(rr[[26]]),unlist(rr[[27]]),unlist(rr[[28]]),unlist(rr[[29]]),unlist(rr[[30]]),unlist(rr[[31]]),unlist(rr[[32]]),unlist(rr[[33]]))
  auc_te_all=c(unlist(rr[[6]]),unlist(rr[[7]]),unlist(rr[[8]]),unlist(rr[[9]]),unlist(rr[[10]]),unlist(rr[[11]]),unlist(rr[[12]]),unlist(rr[[13]]),unlist(rr[[14]]),unlist(rr[[15]]))
  # cal_train=rr[[34]]
  #  cal_test=rr[[35]]
  #  cal_ex=rr[[36]]
  #  return(list(censoring_rate_train,coeff,cindex,length_coef,bs,auc,censoring_rate_test,distance_train,distance_test,lambda,npasses,devratio,cindex_train,bs_train,cindex_ex,bs_ex,auc_ex,auc_ex_all,auc_te_all,cal_train,cal_test,cal_ex))
  return(list(censoring_rate_train,coeff,cindex,length_coef,bs,auc,censoring_rate_test,distance_train,distance_test,lambda,npasses,devratio,cindex_train,bs_train,cindex_ex,bs_ex,auc_ex,auc_ex_all,auc_te_all))
},otherwise = NA)

#write into one function to run in parallell
stability_selection_fun7=purrr::possibly(function(i,data,ex_validation,timess,prob){
  set.seed(i)
  updated_data=bootstrap_fun3(data,i,desired_proportion_A=prob)
  censoring_rate_train=table(updated_data[[1]]$status)[1]/(table(updated_data[[1]]$status)[1]+table(updated_data[[1]]$status)[2])
  censoring_rate_test=table(updated_data[[2]]$status)[1]/(table(updated_data[[2]]$status)[1]+table(updated_data[[2]]$status)[2])
  updated_data_new1=updated_data[[1]]
  updated_data_new1=updated_data_new1[,-which(colnames(updated_data_new1)=="control_var")]
  updated_data_new2=updated_data[[2]]
  updated_data_new2=updated_data_new2[,-which(colnames(updated_data_new2)=="control_var")] 
  updated_data=list(updated_data_new1,updated_data_new2)
  rr=lasso_feature_selection_withoutcal(updated_data,ex_validation,timess)
  coeff=rr[[1]]
  cindex=rr[[2]]
  length_coef=rr[[3]]
  bs=rr[[4]]
  auc=rr[[5]]
  distance_train=ks.test(updated_data[[1]]$time,data$time)
  distance_test=ks.test(updated_data[[2]]$time,data$time)
  lambda=rr[[16]]
  npasses=rr[[17]]
  devratio=rr[[18]]
  cindex_train=rr[[19]]
  bs_train=rr[[20]]
  cindex_ex=rr[[21]]
  bs_ex=rr[[22]]
  auc_ex=rr[[23]]
  auc_ex_all=c(unlist(rr[[24]]),unlist(rr[[25]]),unlist(rr[[26]]),unlist(rr[[27]]),unlist(rr[[28]]),unlist(rr[[29]]),unlist(rr[[30]]),unlist(rr[[31]]),unlist(rr[[32]]),unlist(rr[[33]]))
  auc_te_all=c(unlist(rr[[6]]),unlist(rr[[7]]),unlist(rr[[8]]),unlist(rr[[9]]),unlist(rr[[10]]),unlist(rr[[11]]),unlist(rr[[12]]),unlist(rr[[13]]),unlist(rr[[14]]),unlist(rr[[15]]))
  # cal_train=rr[[34]]
  # cal_test=rr[[35]]
  # cal_ex=rr[[36]]
  # return(list(censoring_rate_train,coeff,cindex,length_coef,bs,auc,censoring_rate_test,distance_train,distance_test,lambda,npasses,devratio,cindex_train,bs_train,cindex_ex,bs_ex,auc_ex,auc_ex_all,auc_te_all,cal_train,cal_test,cal_ex))
  return(list(censoring_rate_train,coeff,cindex,length_coef,bs,auc,censoring_rate_test,distance_train,distance_test,lambda,npasses,devratio,cindex_train,bs_train,cindex_ex,bs_ex,auc_ex,auc_ex_all,auc_te_all))
},otherwise = NA)


# run results random sapling without replacement
#CESC
current_data=read.csv("CESC.csv")
current_data2=current_data
colnames(current_data2)[(dim(current_data2)[2]-2):dim(current_data2)[2]]=c("control_var","status","time")
current_data2$time=as.numeric(current_data2$time)
current_data2=current_data2[current_data2$time>0,]
current_data2=current_data2[!is.na(current_data2$time),]
dim(current_data)
dim(current_data2)
table(current_data2$control_var)
current_data2$control_var=ifelse(current_data2$control_var=="FEMALE",0,1) #female is group0 as the proportion control
set.seed(230720)
rownum=sample(nrow(current_data2),nrow(current_data2)*0.2,replace = FALSE)
ex_validation=current_data2[rownum,]
current_data2=current_data2[-rownum,]
dim(ex_validation)
dim(current_data2)
table(current_data2$status)
summary(current_data2$time)
#define the data with control_var to run the with control_var method
current_data3=current_data2
external_validation3=ex_validation
dim(current_data3)
dim(external_validation3)
#get datasets do not include gender to run other methods
current_data2=current_data2[,-which(colnames(current_data2)=="control_var")]
ex_validation=ex_validation[,-which(colnames(ex_validation)=="control_var")]
dim(ex_validation)
dim(current_data2)

timess=seq(as.numeric(summary(current_data2$time)[2]),as.numeric(summary(current_data2$time)[5]),(as.numeric(summary(current_data2$time)[5])-as.numeric(summary(current_data2$time)[2]))/9)


result1=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.1,ex_validation=ex_validation,timess=timess,mc.cores = 15)
result2=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.2,ex_validation=ex_validation,timess=timess,mc.cores = 15)
result3=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.3,ex_validation=ex_validation,timess=timess,mc.cores = 15)
result4=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.4,ex_validation=ex_validation,timess=timess,mc.cores = 15)
result5=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.5,ex_validation=ex_validation,timess=timess,mc.cores = 15)
result6=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.6,ex_validation=ex_validation,timess=timess,mc.cores = 15)
result7=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.7,ex_validation=ex_validation,timess=timess,mc.cores = 15)
result8=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.8,ex_validation=ex_validation,timess=timess,mc.cores = 15)
result9=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.9,ex_validation=ex_validation,timess=timess,mc.cores = 15)
save(result1,result2,result3,result4,result5,result6,result7,result8,result9,file="CESC.RData")

#COAD
current_data=read.csv("COAD.csv")
current_data2=na.omit(current_data)
colnames(current_data2)[(dim(current_data2)[2]-2):dim(current_data2)[2]]=c("control_var","status","time")
current_data2$time=as.numeric(current_data2$time)
current_data2=current_data2[current_data2$time>0,]
current_data2=current_data2[!is.na(current_data2$time),]
dim(current_data)
dim(current_data2)
table(current_data2$control_var)
current_data2$control_var=ifelse(current_data2$control_var=="FEMALE",0,1) #female is group0 as the proportion control
set.seed(230720)
rownum=sample(nrow(current_data2),nrow(current_data2)*0.2,replace = FALSE)
ex_validation=current_data2[rownum,]
current_data2=current_data2[-rownum,]
dim(ex_validation)
dim(current_data2)
table(current_data2$status)
summary(current_data2$time)
#define the data with control_var to run the with control_var method
current_data3=current_data2
external_validation3=ex_validation
dim(current_data3)
dim(external_validation3)
#get datasets do not include gender to run other methods
current_data2=current_data2[,-which(colnames(current_data2)=="control_var")]
ex_validation=ex_validation[,-which(colnames(ex_validation)=="control_var")]
dim(ex_validation)
dim(current_data2)

timess=seq(as.numeric(summary(current_data2$time)[2]),as.numeric(summary(current_data2$time)[5]),(as.numeric(summary(current_data2$time)[5])-as.numeric(summary(current_data2$time)[2]))/9)

stability_selection_fun5(1,data=current_data2,m_percentage=0.1,ex_validation=ex_validation,timess=timess)
result1=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.1,ex_validation=ex_validation,timess=timess,mc.cores = 15)
result2=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.2,ex_validation=ex_validation,timess=timess,mc.cores = 15)
result3=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.3,ex_validation=ex_validation,timess=timess,mc.cores = 15)
result4=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.4,ex_validation=ex_validation,timess=timess,mc.cores = 15)
result5=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.5,ex_validation=ex_validation,timess=timess,mc.cores = 15)
result6=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.6,ex_validation=ex_validation,timess=timess,mc.cores = 15)
result7=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.7,ex_validation=ex_validation,timess=timess,mc.cores = 15)
result8=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.8,ex_validation=ex_validation,timess=timess,mc.cores = 15)
result9=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.9,ex_validation=ex_validation,timess=timess,mc.cores = 15)
save(result1,result2,result3,result4,result5,result6,result7,result8,result9,file="COAD.RData")


# #ESCA:lots of missingness
# current_data=read.csv("ESCA.csv")
# current_data2=na.omit(current_data)
# colnames(current_data2)[(dim(current_data2)[2]-2):dim(current_data2)[2]]=c("control_var","status","time")
# current_data2$time=as.numeric(current_data2$time)
# current_data2=current_data2[current_data2$time>0,]
# current_data2=current_data2[!is.na(current_data2$time),]
# dim(current_data)
# dim(current_data2)
# table(current_data2$control_var)
# current_data2$control_var=ifelse(current_data2$control_var=="FEMALE",0,1) #female is group0 as the proportion control
# set.seed(230720)
# rownum=sample(nrow(current_data2),nrow(current_data2)*0.2,replace = FALSE)
# ex_validation=current_data2[rownum,]
# current_data2=current_data2[-rownum,]
# dim(ex_validation)
# dim(current_data2)
# table(current_data2$status)
# summary(current_data2$time)
# #define the data with control_var to run the with control_var method
# current_data3=current_data2
# external_validation3=ex_validation
# dim(current_data3)
# dim(external_validation3)
# #get datasets do not include gender to run other methods
# current_data2=current_data2[,-which(colnames(current_data2)=="control_var")]
# ex_validation=ex_validation[,-which(colnames(ex_validation)=="control_var")]
# dim(ex_validation)
# dim(current_data2)
# 
# timess=seq(as.numeric(summary(current_data2$time)[2]),as.numeric(summary(current_data2$time)[5]),(as.numeric(summary(current_data2$time)[5])-as.numeric(summary(current_data2$time)[2]))/9)
# 
# 
# result1=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.1,ex_validation=ex_validation,timess=timess,mc.cores = 15)
# result2=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.2,ex_validation=ex_validation,timess=timess,mc.cores = 15)
# result3=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.3,ex_validation=ex_validation,timess=timess,mc.cores = 15)
# result4=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.4,ex_validation=ex_validation,timess=timess,mc.cores = 15)
# result5=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.5,ex_validation=ex_validation,timess=timess,mc.cores = 15)
# result6=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.6,ex_validation=ex_validation,timess=timess,mc.cores = 15)
# result7=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.7,ex_validation=ex_validation,timess=timess,mc.cores = 15)
# result8=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.8,ex_validation=ex_validation,timess=timess,mc.cores = 15)
# result9=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.9,ex_validation=ex_validation,timess=timess,mc.cores = 15)
# save(result1,result2,result3,result4,result5,result6,result7,result8,result9,file="ESCA.RData")


#GBM
current_data=read.csv("GBM.csv")
current_data2=current_data
colnames(current_data2)[(dim(current_data2)[2]-2):dim(current_data2)[2]]=c("control_var","status","time")
current_data2$time=as.numeric(current_data2$time)
current_data2=current_data2[current_data2$time>0,]
current_data2=current_data2[!is.na(current_data2$time),]
dim(current_data)
dim(current_data2)
table(current_data2$control_var)
current_data2$control_var=ifelse(current_data2$control_var=="FEMALE",0,1) #female is group0 as the proportion control
set.seed(230720)
rownum=sample(nrow(current_data2),nrow(current_data2)*0.2,replace = FALSE)
ex_validation=current_data2[rownum,]
current_data2=current_data2[-rownum,]
dim(ex_validation)
dim(current_data2)
table(current_data2$status)
summary(current_data2$time)
#define the data with control_var to run the with control_var method
current_data3=current_data2
external_validation3=ex_validation
dim(current_data3)
dim(external_validation3)
#get datasets do not include gender to run other methods
current_data2=current_data2[,-which(colnames(current_data2)=="control_var")]
ex_validation=ex_validation[,-which(colnames(ex_validation)=="control_var")]
dim(ex_validation)
dim(current_data2)

timess=seq(as.numeric(summary(current_data2$time)[2]),as.numeric(summary(current_data2$time)[5]),(as.numeric(summary(current_data2$time)[5])-as.numeric(summary(current_data2$time)[2]))/9)


result1=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.1,ex_validation=ex_validation,timess=timess,mc.cores = 15)
result2=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.2,ex_validation=ex_validation,timess=timess,mc.cores = 15)
result3=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.3,ex_validation=ex_validation,timess=timess,mc.cores = 15)
result4=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.4,ex_validation=ex_validation,timess=timess,mc.cores = 15)
result5=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.5,ex_validation=ex_validation,timess=timess,mc.cores = 15)
result6=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.6,ex_validation=ex_validation,timess=timess,mc.cores = 15)
result7=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.7,ex_validation=ex_validation,timess=timess,mc.cores = 15)
result8=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.8,ex_validation=ex_validation,timess=timess,mc.cores = 15)
result9=pbmcapply::pbmclapply(1:500, stability_selection_fun5,data=current_data2,m_percentage=0.9,ex_validation=ex_validation,timess=timess,mc.cores = 15)
save(result1,result2,result3,result4,result5,result6,result7,result8,result9,file="GBM.RData")

#generate new brier score and c-index heatmap
load("GBM.RData")

other_table=function(result,data){
  censoring_rate=c()
  cindex=c()
  num_variables=c()
  brier_score=c()
  auc=c()
  distance_stat=c()
  distance_p=c()
  
  censoring_rate_test=c()
  distance_stat_test=c()
  distance_p_test=c()
  
  lambda=c()
  npasses=c()
  devratio=c()
  
  cindex_train=c()
  brier_score_train=c()
  
  cindex_ex=c()
  brier_score_ex=c()
  auc_ex=c()
  
  # auc_all_test=list()
  # auc_all_ex=list()
  for(i in 1:length(result)){
    if(length(result[[i]])>1){
      censoring_rate[i]=result[[i]][[1]]
      cindex[i]=result[[i]][[3]]
      num_variables[i]=result[[i]][[4]]/dim(data)[2]
      brier_score[i]=result[[i]][[5]]
      auc[i]=result[[i]][[6]]
      distance_stat[i]=result[[i]][[8]]$statistic
      distance_p[i]=result[[i]][[8]]$p.value
      
      censoring_rate_test[i]=result[[i]][[7]]
      distance_stat_test[i]=result[[i]][[9]]$statistic
      distance_p_test[i]=result[[i]][[9]]$p.value
      lambda[i]=result[[i]][[10]]
      npasses[i]=result[[i]][[11]]
      devratio[i]=result[[i]][[12]]
      
      cindex_train[i]=result[[i]][[13]]
      brier_score_train[i]=result[[i]][[14]]
      
      cindex_ex[i]=result[[i]][[15]]
      brier_score_ex[i]=result[[i]][[16]]
      auc_ex[i]=result[[i]][[17]]
      
      # auc_all_test[[i]]=result[[i]][[18]]
      # auc_all_ex[[i]]=result[[i]][[19]]
      
    }else{
      censoring_rate[i]=NA
      cindex[i]=NA
      num_variables[i]=NA
      brier_score[i]=NA
      auc[i]=NA
      distance_stat[i]=NA
      distance_p[i]=NA
      censoring_rate_test[i]=NA
      distance_stat_test[i]=NA
      distance_p_test[i]=NA
      lambda[i]=NA
      npasses[i]=NA
      devratio[i]=NA
      cindex_train[i]=NA
      brier_score_train[i]=NA
      cindex_ex[i]=NA
      brier_score_ex[i]=NA
      auc_ex[i]=NA
      
      # auc_all_test[[i]]=NA
      # auc_all_ex[[i]]=NA
    }}
  return(c(mean(censoring_rate,na.rm=TRUE),mean(num_variables,na.rm=TRUE),mean(cindex,na.rm=TRUE),mean(auc,na.rm=TRUE),mean(brier_score,na.rm=TRUE),mean(distance_stat,na.rm=TRUE),mean(distance_p,na.rm=TRUE),mean(censoring_rate_test,na.rm=TRUE),mean(distance_stat_test,na.rm=TRUE),mean(distance_p_test,na.rm=TRUE),mean(lambda,na.rm=TRUE),mean(npasses,na.rm=TRUE),mean(devratio,na.rm=TRUE),mean(cindex_train,na.rm=TRUE),mean(brier_score_train,na.rm=TRUE),mean(cindex_ex,na.rm=TRUE),mean(brier_score_ex,na.rm=TRUE),mean(auc_ex,na.rm=TRUE)))
}


load("GBM.RData")
heatmap_data=data.frame(nrow=9*18,ncol=18)
heatmap_data[1:162,1]=rep(c("m0.1","m0.2","m0.3","m0.4","m0.5","m0.6","m0.7","m0.8","m0.9"),each=18)
heatmap_data[1:162,2]=rep(c("censoring_rate","num_variables","cindex","auc","brier_score","distance_stat","distance_p","censoring_rate_test","distance_stat_test","distance_p_test","lambda","npasses","devratio","cindex_train","brier_score_train","cindex_ex","brier_score_ex","auc_ex"),9)
heatmap_data[1:162,3]=c(other_table(result1,current_data2),other_table(result2,current_data2),other_table(result3,current_data2),other_table(result4,current_data2),other_table(result5,current_data2),other_table(result6,current_data2),other_table(result7,current_data2),other_table(result8,current_data2),other_table(result9,current_data2))
colnames(heatmap_data)=c("case","metric","value")
heatmap_data$metric=factor(heatmap_data$metric,levels = c("censoring_rate","num_variables","cindex","auc","brier_score","distance_stat","distance_p","censoring_rate_test","distance_stat_test","distance_p_test","lambda","npasses","devratio","cindex_train","brier_score_train","cindex_ex","brier_score_ex","auc_ex"))
gbm_result=heatmap_data

load("COAD.RData")
heatmap_data=data.frame(nrow=9*18,ncol=18)
heatmap_data[1:162,1]=rep(c("m0.1","m0.2","m0.3","m0.4","m0.5","m0.6","m0.7","m0.8","m0.9"),each=18)
heatmap_data[1:162,2]=rep(c("censoring_rate","num_variables","cindex","auc","brier_score","distance_stat","distance_p","censoring_rate_test","distance_stat_test","distance_p_test","lambda","npasses","devratio","cindex_train","brier_score_train","cindex_ex","brier_score_ex","auc_ex"),9)
heatmap_data[1:162,3]=c(other_table(result1,current_data2),other_table(result2,current_data2),other_table(result3,current_data2),other_table(result4,current_data2),other_table(result5,current_data2),other_table(result6,current_data2),other_table(result7,current_data2),other_table(result8,current_data2),other_table(result9,current_data2))
colnames(heatmap_data)=c("case","metric","value")
heatmap_data$metric=factor(heatmap_data$metric,levels = c("censoring_rate","num_variables","cindex","auc","brier_score","distance_stat","distance_p","censoring_rate_test","distance_stat_test","distance_p_test","lambda","npasses","devratio","cindex_train","brier_score_train","cindex_ex","brier_score_ex","auc_ex"))
coad_result=heatmap_data

load("CESC.RData")
heatmap_data=data.frame(nrow=9*18,ncol=18)
heatmap_data[1:162,1]=rep(c("m0.1","m0.2","m0.3","m0.4","m0.5","m0.6","m0.7","m0.8","m0.9"),each=18)
heatmap_data[1:162,2]=rep(c("censoring_rate","num_variables","cindex","auc","brier_score","distance_stat","distance_p","censoring_rate_test","distance_stat_test","distance_p_test","lambda","npasses","devratio","cindex_train","brier_score_train","cindex_ex","brier_score_ex","auc_ex"),9)
heatmap_data[1:162,3]=c(other_table(result1,current_data2),other_table(result2,current_data2),other_table(result3,current_data2),other_table(result4,current_data2),other_table(result5,current_data2),other_table(result6,current_data2),other_table(result7,current_data2),other_table(result8,current_data2),other_table(result9,current_data2))
colnames(heatmap_data)=c("case","metric","value")
heatmap_data$metric=factor(heatmap_data$metric,levels = c("censoring_rate","num_variables","cindex","auc","brier_score","distance_stat","distance_p","censoring_rate_test","distance_stat_test","distance_p_test","lambda","npasses","devratio","cindex_train","brier_score_train","cindex_ex","brier_score_ex","auc_ex"))
cesc_result=heatmap_data


all_result=rbind.data.frame(gbm_result,coad_result,cesc_result)
all_result$category=c(rep("gbm_result",nrow(gbm_result)),rep("coad_result",nrow(coad_result)),rep("cesc_result",nrow(cesc_result)))

heatmap_data=all_result[all_result$case%in%c("m0.2","m0.3","m0.4","m0.5","m0.6","m0.7","m0.8")&all_result$metric%in%c("cindex_train","brier_score_train","censoring_rate","distance_stat"),]
gg_heatmap_data=heatmap_data
#gg_heatmap_data[gg_heatmap_data$metric=="brier_score","value"]=gg_heatmap_data[gg_heatmap_data$metric=="brier_score","value"]/max(gg_heatmap_data[gg_heatmap_data$metric=="brier_score","value"])
gg_heatmap_data$name=paste(gg_heatmap_data$category,gg_heatmap_data$case,sep = "_")
ggheatmap <-ggplot(gg_heatmap_data, aes(name, metric, fill= value)) + 
  geom_tile(color = "grey")+scale_fill_gradientn(limits = c(0,1),colours = c("#2166AC", "#67A9CF" ,"#D1E5F0", "#FFFFFF","#FDDBC7", "#EF8A62", "#B2182B"))+ theme(aspect.ratio = 1, text = element_text(size = 15), legend.position = "bottom") + labs(y= 'Metrics', x = "Cases", fill = 'Values')+theme_bw()+theme(axis.text.x = element_text(color = "black", size = 15,angle = 90),axis.text.y = element_text(color = "black", size = 16))+
  geom_text(aes(label = round(value,3)), color = "black", size = 4)
ggheatmap
ggsave("all_variability_heatmap_train_new.pdf",ggheatmap,width=15,height=10)


heatmap_data=all_result[all_result$case%in%c("m0.2","m0.3","m0.4","m0.5","m0.6","m0.7","m0.8")&all_result$metric%in%c("cindex","brier_score","censoring_rate_test","distance_stat_test"),]
gg_heatmap_data=heatmap_data
#gg_heatmap_data[gg_heatmap_data$metric=="brier_score","value"]=gg_heatmap_data[gg_heatmap_data$metric=="brier_score","value"]/max(gg_heatmap_data[gg_heatmap_data$metric=="brier_score","value"])
gg_heatmap_data$name=paste(gg_heatmap_data$category,gg_heatmap_data$case,sep = "_")
gg_heatmap_data$metric=factor(gg_heatmap_data$metric,levels=rev(c("brier_score","cindex","distance_stat_test","censoring_rate_test")))
ggheatmap <-ggplot(gg_heatmap_data, aes(name, metric, fill= value)) + 
  geom_tile(color = "grey")+scale_fill_gradientn(limits = c(0,1),colours = c("#2166AC", "#67A9CF" ,"#D1E5F0", "#FFFFFF","#FDDBC7", "#EF8A62", "#B2182B"))+ theme(aspect.ratio = 1, text = element_text(size = 15), legend.position = "bottom") + labs(y= 'Metrics', x = "Cases", fill = 'Values')+theme_bw()+theme(axis.text.x = element_text(color = "black", size = 15,angle = 90),axis.text.y = element_text(color = "black", size = 16))+
  geom_text(aes(label = round(value,3)), color = "black", size = 4)
ggheatmap
ggsave("all_variability_heatmap_test_new.pdf",ggheatmap,width=15,height=10)


