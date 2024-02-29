#get boostrapped samples
bootstrap_fun=function(data,i){
  set.seed(i)
  sample_id=sample(1:nrow(data),nrow(data),replace=TRUE)
  bootstrap_sample=data[sample_id,]
  testing=data[-sample_id,]
  return(list(bootstrap_sample,testing))
}

#sample with different censoring rate
bootstrap_fun2=function(data,i,desired_proportion_A){
  set.seed(i)
  # Desired proportion for group A
  #desired_proportion_A <- 0.7
  # Identify the indices of group A
  group_A_indices <- which(data$status == 0)
  # Calculate the number of samples needed for group A
  num_samples_A <- round(desired_proportion_A * length(data$status))
  # Sample from group A
  sampled_indices_A <- sample(group_A_indices, size = num_samples_A, replace = TRUE)
  # Sample the remaining from the other groups
  sampled_indices_others <- sample(setdiff(1:nrow(data),group_A_indices),
                                   size = length(data$status) - num_samples_A, replace = TRUE)
  # Combine the sampled indices
  sampled_indices <- c(sampled_indices_A, sampled_indices_others)
  # Sample the data using the sampled indices
  data[sampled_indices_A,"status"]
  data[sampled_indices_others,"status"]
  sampled_data <- data[sampled_indices, ]
  table(sampled_data$status)
  testing=data[-sampled_indices,]
  return(list(sampled_data,testing))
}

#sample with different another variable: control the group0 as the reference group
bootstrap_fun3=function(data,i,desired_proportion_A){
  #notice here the control_var needs to be coded as 0 and 1
  set.seed(i)
  # Desired proportion for group A
  #desired_proportion_A <- 0.7
  # Identify the indices of group A
  group_A_indices <- which(data$control_var == 0)
  # Calculate the number of samples needed for group A
  num_samples_A <- round(desired_proportion_A * length(data$control_var))
  # Sample from group A
  sampled_indices_A <- sample(group_A_indices, size = num_samples_A, replace = TRUE)
  # Sample the remaining from the other groups
  sampled_indices_others <- sample(setdiff(1:nrow(data),group_A_indices),
                                   size = length(data$control_var) - num_samples_A, replace = TRUE)
  # Combine the sampled indices
  sampled_indices <- c(sampled_indices_A, sampled_indices_others)
  # Sample the data using the sampled indices
  data[sampled_indices_A,"control_var"]
  data[sampled_indices_others,"control_var"]
  sampled_data <- data[sampled_indices, ]
  table(sampled_data$control_var)
  testing=data[-sampled_indices,]
  return(list(sampled_data,testing))
}

moutofn_bootstrap_fun=function(data,i,m_percentage){
  set.seed(i)
  m_number=round(m_percentage*nrow(data),0)
  sampleid=sample(1:nrow(data),m_number,replace=FALSE)
  bootstrap_sample=data[sampleid,]#without replacement
  testing=data[-sampleid,]
  return(list(bootstrap_sample,testing))
}

moutofn_bootstrap_fun_wr=function(data,i,m_percentage){
  set.seed(i)
  m_number=round(m_percentage*nrow(data),0)
  sampleid=sample(1:nrow(data),m_number,replace=TRUE)
  bootstrap_sample=data[sampleid,]#with replacement
  testing=data[-sampleid,]
  return(list(bootstrap_sample,testing))
}

#add in cv and other methods
cv_fun=function(data,i,cvK){
  set.seed(i)
  return_data=list()
  for(j in 1:cvK){
  cvSets = cvTools::cvFolds(nrow(data), cvK) 
  test_id = cvSets$subsets[cvSets$which == j]
  test = data[test_id, ]
  train = data[-test_id, ]
  return_data[[j]]=list(train,test)}
  return(return_data)
}
#moutofn bootstrap with control on censoring rate and distance statistics
moutofn_bootstrap_fun2=function(data,i,m_percentage){
  set.seed(i)
  m_number=round(m_percentage*nrow(data),0)
  sampleid=sample(1:nrow(data),m_number,replace=FALSE)
  bootstrap_sample=data[sampleid,]#without replacement
  testing=data[-sampleid,]
  censoring_rate=table(bootstrap_sample$status)[1]/(table(bootstrap_sample$status)[1]+table(bootstrap_sample$status)[2])
  distance=ks.test(bootstrap_sample$time,data$time)
  if (censoring_rate<0.2|censoring_rate>0.8|distance$p.value<0.05){return(list(NA,NA))}else{
  return(list(bootstrap_sample,testing))}
}

moutofn_bootstrap_fun_wr2=function(data,i,m_percentage){
  set.seed(i)
  m_number=round(m_percentage*nrow(data),0)
  sampleid=sample(1:nrow(data),m_number,replace=TRUE)
  bootstrap_sample=data[sampleid,]#with replacement
  testing=data[-sampleid,]
  censoring_rate=table(bootstrap_sample$status)[1]/(table(bootstrap_sample$status)[1]+table(bootstrap_sample$status)[2])
  distance=ks.test(bootstrap_sample$time,data$time)
  if (censoring_rate<0.2|censoring_rate>0.8|distance$p.value<0.05){return(list(NA,NA))}else{
    return(list(bootstrap_sample,testing))}
}

moutofn_bootstrap_fun3=function(data,i,m_percentage){
  set.seed(i)
  m_number=round(m_percentage*nrow(data),0)
  sampleid=sample(1:nrow(data),m_number,replace=FALSE)
  bootstrap_sample=data[sampleid,]#without replacement
  testing=data[-sampleid,]
  censoring_rate=table(bootstrap_sample$status)[1]/(table(bootstrap_sample$status)[1]+table(bootstrap_sample$status)[2])
  distance=ks.test(bootstrap_sample$time,data$time)
  if (censoring_rate<0.3|censoring_rate>0.7|distance$p.value<0.05){return(list(NA,NA))}else{
    return(list(bootstrap_sample,testing))}
}

moutofn_bootstrap_fun_wr3=function(data,i,m_percentage){
  set.seed(i)
  m_number=round(m_percentage*nrow(data),0)
  sampleid=sample(1:nrow(data),m_number,replace=TRUE)
  bootstrap_sample=data[sampleid,]#with replacement
  testing=data[-sampleid,]
  censoring_rate=table(bootstrap_sample$status)[1]/(table(bootstrap_sample$status)[1]+table(bootstrap_sample$status)[2])
  distance=ks.test(bootstrap_sample$time,data$time)
  if (censoring_rate<0.3|censoring_rate>0.7|distance$p.value<0.05){return(list(NA,NA))}else{
    return(list(bootstrap_sample,testing))}
}

moutofn_bootstrap_fun4=function(data,i,m_percentage){
  set.seed(i)
  m_number=round(m_percentage*nrow(data),0)
  sampleid=sample(1:nrow(data),m_number,replace=FALSE)
  bootstrap_sample=data[sampleid,]#without replacement
  testing=data[-sampleid,]
  censoring_rate=table(bootstrap_sample$status)[1]/(table(bootstrap_sample$status)[1]+table(bootstrap_sample$status)[2])
  distance=ks.test(bootstrap_sample$time,data$time)
  if (censoring_rate<(table(current_data2$status)[1]/dim(current_data2)[1]-0.1)|censoring_rate>(table(current_data2$status)[1]/dim(current_data2)[1]+0.1)|distance$p.value<0.05){return(list(NA,NA))}else{
    return(list(bootstrap_sample,testing))}
}

moutofn_bootstrap_fun_wr4=function(data,i,m_percentage){
  set.seed(i)
  m_number=round(m_percentage*nrow(data),0)
  sampleid=sample(1:nrow(data),m_number,replace=TRUE)
  bootstrap_sample=data[sampleid,]#with replacement
  testing=data[-sampleid,]
  censoring_rate=table(bootstrap_sample$status)[1]/(table(bootstrap_sample$status)[1]+table(bootstrap_sample$status)[2])
  distance=ks.test(bootstrap_sample$time,data$time)
  if (censoring_rate<(table(current_data2$status)[1]/dim(current_data2)[1]-0.1)|censoring_rate>(table(current_data2$status)[1]/dim(current_data2)[1]+0.1)|distance$p.value<0.05){return(list(NA,NA))}else{
    return(list(bootstrap_sample,testing))}
}
########################################
#########################################
######################################

lasso_feature_selection_lessgroup=function(data,ex_validation,timess){
  test = data[[2]]
  train =data[[1]]
  tr_matrix=data.matrix(train[,!names(train)%in% c("time","status")])
  te_matrix=data.matrix(test[,!names(test)%in% c("time","status")])
  tr_st=train$time
  tr_y=train$status
  te_st=test$time
  te_y=test$status
  
  
  # suppressMessages(library("doParallel"))
  # registerDoParallel(detectCores())
  fit <- fit_lasso(tr_matrix, Surv(tr_st, tr_y), nfolds = 5, rule = "lambda.min")
  #store the model info
  lambda=fit$model$lambda
  npasses=fit$model$npasses
  devratio=fit$model$dev.ratio
  #names(fit)
  coeff=infer_variable_type(fit, tr_matrix)$name
  if(length(coeff)==0){coeff=NA}
  tr_times=train$time
  te_times=test$time
  pred_te_matrix=predict(fit,tr_matrix, Surv(tr_st, tr_y), newx = te_matrix,pred.at=te_times)
  pred_tr_matrix=predict(fit,tr_matrix, Surv(tr_st, tr_y), newx = tr_matrix,pred.at=tr_times)
  pred_te=c()
  for(i in 1:nrow(pred_te_matrix)){
    pred_te[i]=pred_te_matrix[i,][which(colnames(pred_te_matrix)==te_times[i])[1]]
  }
  pred_tr=c()
  for(i in 1:nrow(pred_tr_matrix)){
    pred_tr[i]=pred_tr_matrix[i,][which(colnames(pred_tr_matrix)==tr_times[i])[1]]
  }
  
  harrelC1 <- rcorr.cens(-pred_tr,with(train,Surv(time,status)))
  cindex_train<-harrelC1["C Index"]
  harrelC1 <- rcorr.cens(-pred_te,with(test,Surv(time,status)))
  cindex_test<-harrelC1["C Index"]
  # #calculate brier score:currently cannot be calculated if only have training data
  times=median(train$time)
  Surv.rsp <- Surv(train$time,train$status)
  Surv.rsp.new <- Surv(test$time,test$status)
  lp<-pred_tr
  lp.new=pred_te
  bs_test <- predErr(Surv.rsp,Surv.rsp.new,lp,lp.new, time = times)$error
  #maybe this is ok for the training data only (naive method)
  calculate_brier_score <- function(predicted_probs, event_status) {
    n <- length(predicted_probs)
    brier_score <- sum((predicted_probs - event_status)^2) / n
    return(brier_score)
  }
  bs_train=calculate_brier_score(lp,train$status)
  #calculate time-dependent AUC
  times=timess
  auc.uno=AUC.uno(Surv.rsp, Surv.rsp.new, lp.new, times)
  a1=auc.uno$auc[1]
  a2=auc.uno$auc[2]
  a3=auc.uno$auc[3]
  a4=auc.uno$auc[4]
  a5=auc.uno$auc[5]
  a6=auc.uno$auc[6]
  a7=auc.uno$auc[7]
  a8=auc.uno$auc[8]
  a9=auc.uno$auc[9]
  a10=auc.uno$auc[10]
  auc=auc.uno$iauc
  
  # Calculate the calibration slope:from packages rms, pec, riskregression, they all have the calibration function, however, calibration plot is dependent of the model fit, so it cannot be obtained with all models
  #calculate the calibration plot within this package
  cal_train=calibrate(tr_matrix,tr_st,tr_y,model.type = "lasso",pred.at = median(train$time),alpha=1,lambda = fit$lambda,method = "fitting",ngroup = 2,seed = 1010)
  cal_test=calibrate_external(object=fit,tr_matrix,tr_st,tr_y,x_new=te_matrix,time_new=te_st, event_new=te_y,pred.at = median(train$time),ngroup = 2)
  
  #testing results on external validation datasets
  ex_matrix=data.matrix(ex_validation[,!names(ex_validation)%in% c("time","status")])
  ex_st=ex_validation$time
  ex_y=ex_validation$status
  ex_times=ex_validation$time
  pred_ex_matrix=predict(fit,ex_matrix, Surv(ex_st, ex_y), newx = ex_matrix,pred.at=ex_times)
  pred_ex=c()
  for(i in 1:nrow(pred_ex_matrix)){
    pred_ex[i]=pred_ex_matrix[i,][which(colnames(pred_ex_matrix)==ex_times[i])[1]]
  }
  harrelC1 <- rcorr.cens(-pred_ex,with(ex_validation,Surv(time,status)))
  cindex_ex<-harrelC1["C Index"]
  
  Surv.rsp <- Surv(train$time,train$status)
  Surv.rsp.new <- Surv(ex_validation$time,ex_validation$status)
  lp<-pred_tr
  lp.new=pred_ex
  times=median(train$time)
  bs_ex <- predErr(Surv.rsp,Surv.rsp.new,lp,lp.new, time = times)$error
  #calculate time-dependent AUC
  times=timess
  auc.uno=AUC.uno(Surv.rsp, Surv.rsp.new, lp.new, times)
  a1_ex=auc.uno$auc[1]
  a2_ex=auc.uno$auc[2]
  a3_ex=auc.uno$auc[3]
  a4_ex=auc.uno$auc[4]
  a5_ex=auc.uno$auc[5]
  a6_ex=auc.uno$auc[6]
  a7_ex=auc.uno$auc[7]
  a8_ex=auc.uno$auc[8]
  a9_ex=auc.uno$auc[9]
  a10_ex=auc.uno$auc[10]
  auc_ex=auc.uno$iauc
  #calculate calibration curve
  cal_ex=calibrate_external(object=fit,tr_matrix,tr_st,tr_y,x_new=ex_matrix,time_new=ex_st, event_new=ex_y,pred.at = median(train$time),ngroup = 2)
  
  return(list(coeff,cindex_test,length(coeff),bs_test,auc,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,lambda,npasses,devratio,cindex_train,bs_train,cindex_ex,bs_ex,auc_ex,a1_ex,a2_ex,a3_ex,a4_ex,a5_ex,a6_ex,a7_ex,a8_ex,a9_ex,a10_ex,cal_train,cal_test,cal_ex))
}

lasso_feature_selection=function(data,ex_validation,timess){
  test = data[[2]]
  train =data[[1]]
  tr_matrix=data.matrix(train[,!names(train)%in% c("time","status")])
  te_matrix=data.matrix(test[,!names(test)%in% c("time","status")])
  tr_st=train$time
  tr_y=train$status
  te_st=test$time
  te_y=test$status
  
  
  # suppressMessages(library("doParallel"))
  # registerDoParallel(detectCores())
  fit <- fit_lasso(tr_matrix, Surv(tr_st, tr_y), nfolds = 5, rule = "lambda.min")
  #store the model info
  lambda=fit$model$lambda
  npasses=fit$model$npasses
  devratio=fit$model$dev.ratio
  #names(fit)
  coeff=infer_variable_type(fit, tr_matrix)$name
  if(length(coeff)==0){coeff=NA}
  tr_times=train$time
  te_times=test$time
  pred_te_matrix=predict(fit,tr_matrix, Surv(tr_st, tr_y), newx = te_matrix,pred.at=te_times)
  pred_tr_matrix=predict(fit,tr_matrix, Surv(tr_st, tr_y), newx = tr_matrix,pred.at=tr_times)
  pred_te=c()
  for(i in 1:nrow(pred_te_matrix)){
    pred_te[i]=pred_te_matrix[i,][which(colnames(pred_te_matrix)==te_times[i])[1]]
  }
  pred_tr=c()
  for(i in 1:nrow(pred_tr_matrix)){
    pred_tr[i]=pred_tr_matrix[i,][which(colnames(pred_tr_matrix)==tr_times[i])[1]]
  }
  
  harrelC1 <- rcorr.cens(-pred_tr,with(train,Surv(time,status)))
  cindex_train<-harrelC1["C Index"]
  harrelC1 <- rcorr.cens(-pred_te,with(test,Surv(time,status)))
  cindex_test<-harrelC1["C Index"]
  # #calculate brier score:currently cannot be calculated if only have training data
  times=median(train$time)
  Surv.rsp <- Surv(train$time,train$status)
  Surv.rsp.new <- Surv(test$time,test$status)
  lp<-pred_tr
  lp.new=pred_te
  bs_test <- predErr(Surv.rsp,Surv.rsp.new,lp,lp.new, time = times)$error
  #maybe this is ok for the training data only (naive method)
  calculate_brier_score <- function(predicted_probs, event_status) {
    n <- length(predicted_probs)
    brier_score <- sum((predicted_probs - event_status)^2) / n
    return(brier_score)
  }
  bs_train=calculate_brier_score(lp,train$status)
  #calculate time-dependent AUC
  times=timess
  auc.uno=AUC.uno(Surv.rsp, Surv.rsp.new, lp.new, times)
  a1=auc.uno$auc[1]
  a2=auc.uno$auc[2]
  a3=auc.uno$auc[3]
  a4=auc.uno$auc[4]
  a5=auc.uno$auc[5]
  a6=auc.uno$auc[6]
  a7=auc.uno$auc[7]
  a8=auc.uno$auc[8]
  a9=auc.uno$auc[9]
  a10=auc.uno$auc[10]
  auc=auc.uno$iauc
  
  # Calculate the calibration slope:from packages rms, pec, riskregression, they all have the calibration function, however, calibration plot is dependent of the model fit, so it cannot be obtained with all models
  #calculate the calibration plot within this package
  cal_train=calibrate(tr_matrix,tr_st,tr_y,model.type = "lasso",pred.at = median(train$time),alpha=1,lambda = fit$lambda,method = "fitting",ngroup = 5,seed = 1010)
  cal_test=calibrate_external(object=fit,tr_matrix,tr_st,tr_y,x_new=te_matrix,time_new=te_st, event_new=te_y,pred.at = median(train$time),ngroup = 5)
  
  #testing results on external validation datasets
  ex_matrix=data.matrix(ex_validation[,!names(ex_validation)%in% c("time","status")])
  ex_st=ex_validation$time
  ex_y=ex_validation$status
  ex_times=ex_validation$time
  pred_ex_matrix=predict(fit,ex_matrix, Surv(ex_st, ex_y), newx = ex_matrix,pred.at=ex_times)
  pred_ex=c()
  for(i in 1:nrow(pred_ex_matrix)){
    pred_ex[i]=pred_ex_matrix[i,][which(colnames(pred_ex_matrix)==ex_times[i])[1]]
  }
  harrelC1 <- rcorr.cens(-pred_ex,with(ex_validation,Surv(time,status)))
  cindex_ex<-harrelC1["C Index"]
  
  Surv.rsp <- Surv(train$time,train$status)
  Surv.rsp.new <- Surv(ex_validation$time,ex_validation$status)
  lp<-pred_tr
  lp.new=pred_ex
  times=median(train$time)
  bs_ex <- predErr(Surv.rsp,Surv.rsp.new,lp,lp.new, time = times)$error
  #calculate time-dependent AUC
  times=timess
  auc.uno=AUC.uno(Surv.rsp, Surv.rsp.new, lp.new, times)
  a1_ex=auc.uno$auc[1]
  a2_ex=auc.uno$auc[2]
  a3_ex=auc.uno$auc[3]
  a4_ex=auc.uno$auc[4]
  a5_ex=auc.uno$auc[5]
  a6_ex=auc.uno$auc[6]
  a7_ex=auc.uno$auc[7]
  a8_ex=auc.uno$auc[8]
  a9_ex=auc.uno$auc[9]
  a10_ex=auc.uno$auc[10]
  auc_ex=auc.uno$iauc
  #calculate calibration curve
  cal_ex=calibrate_external(object=fit,tr_matrix,tr_st,tr_y,x_new=ex_matrix,time_new=ex_st, event_new=ex_y,pred.at = median(train$time),ngroup = 5)
  
  return(list(coeff,cindex_test,length(coeff),bs_test,auc,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,lambda,npasses,devratio,cindex_train,bs_train,cindex_ex,bs_ex,auc_ex,a1_ex,a2_ex,a3_ex,a4_ex,a5_ex,a6_ex,a7_ex,a8_ex,a9_ex,a10_ex,cal_train,cal_test,cal_ex))
}

lasso_feature_selection_withoutcal=function(data,ex_validation,timess){
  test = data[[2]]
  train =data[[1]]
  tr_matrix=data.matrix(train[,!names(train)%in% c("time","status")])
  te_matrix=data.matrix(test[,!names(test)%in% c("time","status")])
  tr_st=train$time
  tr_y=train$status
  te_st=test$time
  te_y=test$status
  
  
  # suppressMessages(library("doParallel"))
  # registerDoParallel(detectCores())
  fit <- fit_lasso(tr_matrix, Surv(tr_st, tr_y), nfolds = 5, rule = "lambda.min")
  #store the model info
  lambda=fit$model$lambda
  npasses=fit$model$npasses
  devratio=fit$model$dev.ratio
  #names(fit)
  coeff=infer_variable_type(fit, tr_matrix)$name
  if(length(coeff)==0){coeff=NA}
  tr_times=train$time
  te_times=test$time
  pred_te_matrix=predict(fit,tr_matrix, Surv(tr_st, tr_y), newx = te_matrix,pred.at=te_times)
  pred_tr_matrix=predict(fit,tr_matrix, Surv(tr_st, tr_y), newx = tr_matrix,pred.at=tr_times)
  pred_te=c()
  for(i in 1:nrow(pred_te_matrix)){
    pred_te[i]=pred_te_matrix[i,][which(colnames(pred_te_matrix)==te_times[i])[1]]
  }
  pred_tr=c()
  for(i in 1:nrow(pred_tr_matrix)){
    pred_tr[i]=pred_tr_matrix[i,][which(colnames(pred_tr_matrix)==tr_times[i])[1]]
  }
  
  harrelC1 <- rcorr.cens(-pred_tr,with(train,Surv(time,status)))
  cindex_train<-harrelC1["C Index"]
  harrelC1 <- rcorr.cens(-pred_te,with(test,Surv(time,status)))
  cindex_test<-harrelC1["C Index"]
  # #calculate brier score:currently cannot be calculated if only have training data
  times=median(train$time)
  Surv.rsp <- Surv(train$time,train$status)
  Surv.rsp.new <- Surv(test$time,test$status)
  lp<-pred_tr
  lp.new=pred_te
  bs_test <- predErr(Surv.rsp,Surv.rsp.new,lp,lp.new, time = times)$error
  #maybe this is ok for the training data only (naive method)
  calculate_brier_score <- function(predicted_probs, event_status) {
    n <- length(predicted_probs)
    brier_score <- sum((predicted_probs - event_status)^2) / n
    return(brier_score)
  }
  bs_train=calculate_brier_score(lp,train$status)
  #calculate time-dependent AUC
  times=timess
  auc.uno=AUC.uno(Surv.rsp, Surv.rsp.new, lp.new, times)
  a1=auc.uno$auc[1]
  a2=auc.uno$auc[2]
  a3=auc.uno$auc[3]
  a4=auc.uno$auc[4]
  a5=auc.uno$auc[5]
  a6=auc.uno$auc[6]
  a7=auc.uno$auc[7]
  a8=auc.uno$auc[8]
  a9=auc.uno$auc[9]
  a10=auc.uno$auc[10]
  auc=auc.uno$iauc
  
  # # Calculate the calibration slope:from packages rms, pec, riskregression, they all have the calibration function, however, calibration plot is dependent of the model fit, so it cannot be obtained with all models
  # #calculate the calibration plot within this package
  # cal_train=calibrate(tr_matrix,tr_st,tr_y,model.type = "lasso",pred.at = median(train$time),alpha=1,lambda = fit$lambda,method = "fitting",ngroup = 5,seed = 1010)
  # cal_test=calibrate_external(object=fit,tr_matrix,tr_st,tr_y,x_new=te_matrix,time_new=te_st, event_new=te_y,pred.at = median(train$time),ngroup = 5)
  # 
  #testing results on external validation datasets
  ex_matrix=data.matrix(ex_validation[,!names(ex_validation)%in% c("time","status")])
  ex_st=ex_validation$time
  ex_y=ex_validation$status
  ex_times=ex_validation$time
  pred_ex_matrix=predict(fit,ex_matrix, Surv(ex_st, ex_y), newx = ex_matrix,pred.at=ex_times)
  pred_ex=c()
  for(i in 1:nrow(pred_ex_matrix)){
    pred_ex[i]=pred_ex_matrix[i,][which(colnames(pred_ex_matrix)==ex_times[i])[1]]
  }
  harrelC1 <- rcorr.cens(-pred_ex,with(ex_validation,Surv(time,status)))
  cindex_ex<-harrelC1["C Index"]
  
  Surv.rsp <- Surv(train$time,train$status)
  Surv.rsp.new <- Surv(ex_validation$time,ex_validation$status)
  lp<-pred_tr
  lp.new=pred_ex
  times=median(train$time)
  bs_ex <- predErr(Surv.rsp,Surv.rsp.new,lp,lp.new, time = times)$error
  #calculate time-dependent AUC
  times=timess
  auc.uno=AUC.uno(Surv.rsp, Surv.rsp.new, lp.new, times)
  a1_ex=auc.uno$auc[1]
  a2_ex=auc.uno$auc[2]
  a3_ex=auc.uno$auc[3]
  a4_ex=auc.uno$auc[4]
  a5_ex=auc.uno$auc[5]
  a6_ex=auc.uno$auc[6]
  a7_ex=auc.uno$auc[7]
  a8_ex=auc.uno$auc[8]
  a9_ex=auc.uno$auc[9]
  a10_ex=auc.uno$auc[10]
  auc_ex=auc.uno$iauc
  #calculate calibration curve
  #cal_ex=calibrate_external(object=fit,tr_matrix,tr_st,tr_y,x_new=ex_matrix,time_new=ex_st, event_new=ex_y,pred.at = median(train$time),ngroup = 5)
  
  return(list(coeff,cindex_test,length(coeff),bs_test,auc,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,lambda,npasses,devratio,cindex_train,bs_train,cindex_ex,bs_ex,auc_ex,a1_ex,a2_ex,a3_ex,a4_ex,a5_ex,a6_ex,a7_ex,a8_ex,a9_ex,a10_ex))
}

#timess=seq(as.numeric(summary(current_data2$time)[2]),as.numeric(summary(current_data2$time)[5]),(as.numeric(summary(current_data2$time)[5])-as.numeric(summary(current_data2$time)[2]))/9)

lasso_feature_selection2=function(data,ex_validation,timess){
  test = data[[2]]
  train =data[[1]]
  tr_matrix=data.matrix(train[,!names(train)%in% c("time","status")])
  te_matrix=data.matrix(test[,!names(test)%in% c("time","status")])
  tr_st=train$time
  tr_y=train$status
  te_st=test$time
  te_y=test$status
  
  
  fit0=cv.glmnet(tr_matrix, Surv(tr_st,tr_y), family="cox", standardize = F,alpha=1, nfolds = 5,type.measure = "C")
  fit=glmnet(tr_matrix, Surv(tr_st,tr_y), family="cox", standardize = F,alpha=1,lambda = fit0$lambda.min,type.measure = "C")
  #store the model info
  lambda=fit$lambda
  npasses=fit$npasses
  devratio=fit$dev.ratio
  #names(fit)
  coeff=rownames(coef(fit, s = 'lambda.min'))[coef(fit, s = 'lambda.min')[,1]!= 0]
  if(length(coeff)==0){coeff=NA}
  
  pred_tr<-predict(fit,newx=tr_matrix)
  pred_te=predict(fit,newx=te_matrix)
  
  harrelC1 <- rcorr.cens(-pred_tr,with(train,Surv(time,status)))
  cindex_train<-harrelC1["C Index"]
  harrelC1 <- rcorr.cens(-pred_te,with(test,Surv(time,status)))
  cindex_test<-harrelC1["C Index"]
  # #calculate brier score:currently cannot be calculated if only have training data
  times=median(train$time)
  Surv.rsp <- Surv(train$time,train$status)
  Surv.rsp.new <- Surv(test$time,test$status)
  lp<-pred_tr
  lp.new=pred_te
  bs_test <- predErr(Surv.rsp,Surv.rsp.new,lp,lp.new, time = times)$error
  #maybe this is ok for the training data only (naive method)
  calculate_brier_score <- function(predicted_probs, event_status) {
    n <- length(predicted_probs)
    brier_score <- sum((predicted_probs - event_status)^2) / n
    return(brier_score)
  }
  bs_train=calculate_brier_score(lp,train$status)
  #calculate time-dependent AUC
  times=timess
  auc.uno=AUC.uno(Surv.rsp, Surv.rsp.new, lp.new, times)
  a1=auc.uno$auc[1]
  a2=auc.uno$auc[2]
  a3=auc.uno$auc[3]
  a4=auc.uno$auc[4]
  a5=auc.uno$auc[5]
  a6=auc.uno$auc[6]
  a7=auc.uno$auc[7]
  a8=auc.uno$auc[8]
  a9=auc.uno$auc[9]
  a10=auc.uno$auc[10]
  auc=auc.uno$iauc
  
  # Calculate the calibration slope:from packages rms, pec, riskregression, they all have the calibration function, however, calibration plot is dependent of the model fit, so it cannot be obtained with all models
  #calculate a naive calibration slope for the prediction
  
  #testing results on external validation datasets
  ex_matrix=data.matrix(ex_validation[,!names(ex_validation)%in% c("time","status")])
  ex_st=ex_validation$time
  ex_y=ex_validation$status
  
  pred_ex=predict(fit,newx=ex_matrix)
  harrelC1 <- rcorr.cens(-pred_ex,with(ex_validation,Surv(time,status)))
  cindex_ex<-harrelC1["C Index"]
  
  Surv.rsp <- Surv(train$time,train$status)
  Surv.rsp.new <- Surv(ex_validation$time,ex_validation$status)
  lp<-pred_tr
  lp.new=pred_ex
  times=median(train$time)
  bs_ex <- predErr(Surv.rsp,Surv.rsp.new,lp,lp.new, time = times)$error
  #calculate time-dependent AUC
  times=timess
  auc.uno=AUC.uno(Surv.rsp, Surv.rsp.new, lp.new, times)
  a1_ex=auc.uno$auc[1]
  a2_ex=auc.uno$auc[2]
  a3_ex=auc.uno$auc[3]
  a4_ex=auc.uno$auc[4]
  a5_ex=auc.uno$auc[5]
  a6_ex=auc.uno$auc[6]
  a7_ex=auc.uno$auc[7]
  a8_ex=auc.uno$auc[8]
  a9_ex=auc.uno$auc[9]
  a10_ex=auc.uno$auc[10]
  auc_ex=auc.uno$iauc
  
  return(list(coeff,cindex_test,length(coeff),bs_test,auc,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,lambda,npasses,devratio,cindex_train,bs_train,cindex_ex,bs_ex,auc_ex,a1_ex,a2_ex,a3_ex,a4_ex,a5_ex,a6_ex,a7_ex,a8_ex,a9_ex,a10_ex))
}


enet_feature_selection2=function(data,ex_validation,timess){
  test = data[[2]]
  train =data[[1]]
  tr_matrix=data.matrix(train[,!names(train)%in% c("time","status")])
  te_matrix=data.matrix(test[,!names(test)%in% c("time","status")])
  tr_st=train$time
  tr_y=train$status
  te_st=test$time
  te_y=test$status
  
  
  fit0=cv.glmnet(tr_matrix, Surv(tr_st,tr_y), family="cox", standardize = F,alpha=0.5, nfolds = 5,type.measure = "C")
  fit=glmnet(tr_matrix, Surv(tr_st,tr_y), family="cox", standardize = F,alpha=0.5,lambda = fit0$lambda.min,type.measure = "C")
  #store the model info
  lambda=fit$lambda
  npasses=fit$npasses
  devratio=fit$dev.ratio
  #names(fit)
  coeff=rownames(coef(fit, s = 'lambda.min'))[coef(fit, s = 'lambda.min')[,1]!= 0]
  if(length(coeff)==0){coeff=NA}
  
  pred_tr<-predict(fit,newx=tr_matrix)
  pred_te=predict(fit,newx=te_matrix)
  
  harrelC1 <- rcorr.cens(-pred_tr,with(train,Surv(time,status)))
  cindex_train<-harrelC1["C Index"]
  harrelC1 <- rcorr.cens(-pred_te,with(test,Surv(time,status)))
  cindex_test<-harrelC1["C Index"]
  # #calculate brier score:currently cannot be calculated if only have training data
  times=median(train$time)
  Surv.rsp <- Surv(train$time,train$status)
  Surv.rsp.new <- Surv(test$time,test$status)
  lp<-pred_tr
  lp.new=pred_te
  bs_test <- predErr(Surv.rsp,Surv.rsp.new,lp,lp.new, time = times)$error
  #maybe this is ok for the training data only (naive method)
  calculate_brier_score <- function(predicted_probs, event_status) {
    n <- length(predicted_probs)
    brier_score <- sum((predicted_probs - event_status)^2) / n
    return(brier_score)
  }
  bs_train=calculate_brier_score(lp,train$status)
  #calculate time-dependent AUC
  times=timess
  auc.uno=AUC.uno(Surv.rsp, Surv.rsp.new, lp.new, times)
  a1=auc.uno$auc[1]
  a2=auc.uno$auc[2]
  a3=auc.uno$auc[3]
  a4=auc.uno$auc[4]
  a5=auc.uno$auc[5]
  a6=auc.uno$auc[6]
  a7=auc.uno$auc[7]
  a8=auc.uno$auc[8]
  a9=auc.uno$auc[9]
  a10=auc.uno$auc[10]
  auc=auc.uno$iauc
  
  # Calculate the calibration slope:from packages rms, pec, riskregression, they all have the calibration function, however, calibration plot is dependent of the model fit, so it cannot be obtained with all models
  #calculate a naive calibration slope for the prediction
  
  #testing results on external validation datasets
  ex_matrix=data.matrix(ex_validation[,!names(ex_validation)%in% c("time","status")])
  ex_st=ex_validation$time
  ex_y=ex_validation$status
  
  pred_ex=predict(fit,newx=ex_matrix)
  harrelC1 <- rcorr.cens(-pred_ex,with(ex_validation,Surv(time,status)))
  cindex_ex<-harrelC1["C Index"]
  
  Surv.rsp <- Surv(train$time,train$status)
  Surv.rsp.new <- Surv(ex_validation$time,ex_validation$status)
  lp<-pred_tr
  lp.new=pred_ex
  times=median(train$time)
  bs_ex <- predErr(Surv.rsp,Surv.rsp.new,lp,lp.new, time = times)$error
  #calculate time-dependent AUC
  times=timess
  auc.uno=AUC.uno(Surv.rsp, Surv.rsp.new, lp.new, times)
  a1_ex=auc.uno$auc[1]
  a2_ex=auc.uno$auc[2]
  a3_ex=auc.uno$auc[3]
  a4_ex=auc.uno$auc[4]
  a5_ex=auc.uno$auc[5]
  a6_ex=auc.uno$auc[6]
  a7_ex=auc.uno$auc[7]
  a8_ex=auc.uno$auc[8]
  a9_ex=auc.uno$auc[9]
  a10_ex=auc.uno$auc[10]
  auc_ex=auc.uno$iauc
  
  return(list(coeff,cindex_test,length(coeff),bs_test,auc,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,lambda,npasses,devratio,cindex_train,bs_train,cindex_ex,bs_ex,auc_ex,a1_ex,a2_ex,a3_ex,a4_ex,a5_ex,a6_ex,a7_ex,a8_ex,a9_ex,a10_ex))
}

enet_feature_selection=function(data,ex_validation,timess){
  test = data[[2]]
  train =data[[1]]
  tr_matrix=data.matrix(train[,!names(train)%in% c("time","status")])
  te_matrix=data.matrix(test[,!names(test)%in% c("time","status")])
  tr_st=train$time
  tr_y=train$status
  te_st=test$time
  te_y=test$status
  
  # suppressMessages(library("doParallel"))
  # registerDoParallel(detectCores())
  fit <- fit_enet(tr_matrix, Surv(tr_st, tr_y), nfolds = 5, rule = "lambda.min")
  #store the model info
  lambda=fit$model$lambda
  npasses=fit$model$npasses
  devratio=fit$model$dev.ratio
  alpha=fit$alpha
  #names(fit)
  coeff=infer_variable_type(fit, tr_matrix)$name
  if(length(coeff)==0){coeff=NA}
  tr_times=train$time
  te_times=test$time
  pred_te_matrix=predict(fit,tr_matrix, Surv(tr_st, tr_y), newx = te_matrix,pred.at=te_times)
  pred_tr_matrix=predict(fit,tr_matrix, Surv(tr_st, tr_y), newx = tr_matrix,pred.at=tr_times)
  pred_te=c()
  for(i in 1:nrow(pred_te_matrix)){
    pred_te[i]=pred_te_matrix[i,][which(colnames(pred_te_matrix)==te_times[i])[1]]
  }
  pred_tr=c()
  for(i in 1:nrow(pred_tr_matrix)){
    pred_tr[i]=pred_tr_matrix[i,][which(colnames(pred_tr_matrix)==tr_times[i])[1]]
  }
  
  harrelC1 <- rcorr.cens(-pred_tr,with(train,Surv(time,status)))
  cindex_train<-harrelC1["C Index"]
  harrelC1 <- rcorr.cens(-pred_te,with(test,Surv(time,status)))
  cindex_test<-harrelC1["C Index"]
  # #calculate brier score:currently cannot be calculated if only have training data
  times=median(train$time)
  Surv.rsp <- Surv(train$time,train$status)
  Surv.rsp.new <- Surv(test$time,test$status)
  lp<-pred_tr
  lp.new=pred_te
  bs_test <- predErr(Surv.rsp,Surv.rsp.new,lp,lp.new, time = times)$error
  #maybe this is ok for the training data only (naive method)
  calculate_brier_score <- function(predicted_probs, event_status) {
    n <- length(predicted_probs)
    brier_score <- sum((predicted_probs - event_status)^2) / n
    return(brier_score)
  }
  bs_train=calculate_brier_score(lp,train$status)
  #calculate time-dependent AUC
  times=timess
  auc.uno=AUC.uno(Surv.rsp, Surv.rsp.new, lp.new, times)
  a1=auc.uno$auc[1]
  a2=auc.uno$auc[2]
  a3=auc.uno$auc[3]
  a4=auc.uno$auc[4]
  a5=auc.uno$auc[5]
  a6=auc.uno$auc[6]
  a7=auc.uno$auc[7]
  a8=auc.uno$auc[8]
  a9=auc.uno$auc[9]
  a10=auc.uno$auc[10]
  auc=auc.uno$iauc
  
  # Calculate the calibration slope:from packages rms, pec, riskregression, they all have the calibration function, however, calibration plot is dependent of the model fit, so it cannot be obtained with all models
  #calculate a naive calibration slope for the prediction
  
  #testing results on external validation datasets
  ex_matrix=data.matrix(ex_validation[,!names(ex_validation)%in% c("time","status")])
  ex_st=ex_validation$time
  ex_y=ex_validation$status
  ex_times=ex_validation$time
  pred_ex_matrix=predict(fit,ex_matrix, Surv(ex_st, ex_y), newx = ex_matrix,pred.at=ex_times)
  pred_ex=c()
  for(i in 1:nrow(pred_ex_matrix)){
    pred_ex[i]=pred_ex_matrix[i,][which(colnames(pred_ex_matrix)==ex_times[i])[1]]
  }
  harrelC1 <- rcorr.cens(-pred_ex,with(ex_validation,Surv(time,status)))
  cindex_ex<-harrelC1["C Index"]
  
  Surv.rsp <- Surv(train$time,train$status)
  Surv.rsp.new <- Surv(ex_validation$time,ex_validation$status)
  lp<-pred_tr
  lp.new=pred_ex
  times=median(train$time)
  bs_ex <- predErr(Surv.rsp,Surv.rsp.new,lp,lp.new, time = times)$error
  #calculate time-dependent AUC
  times=timess
  auc.uno=AUC.uno(Surv.rsp, Surv.rsp.new, lp.new, times)
  a1_ex=auc.uno$auc[1]
  a2_ex=auc.uno$auc[2]
  a3_ex=auc.uno$auc[3]
  a4_ex=auc.uno$auc[4]
  a5_ex=auc.uno$auc[5]
  a6_ex=auc.uno$auc[6]
  a7_ex=auc.uno$auc[7]
  a8_ex=auc.uno$auc[8]
  a9_ex=auc.uno$auc[9]
  a10_ex=auc.uno$auc[10]
  auc_ex=auc.uno$iauc
  
  return(list(coeff,cindex_test,length(coeff),bs_test,auc,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,lambda,npasses,devratio,cindex_train,bs_train,cindex_ex,bs_ex,auc_ex,a1_ex,a2_ex,a3_ex,a4_ex,a5_ex,a6_ex,a7_ex,a8_ex,a9_ex,a10_ex,alpha))
}

alasso_feature_selection=function(data,ex_validation,timess){
  test = data[[2]]
  train =data[[1]]
  tr_matrix=data.matrix(train[,!names(train)%in% c("time","status")])
  te_matrix=data.matrix(test[,!names(test)%in% c("time","status")])
  tr_st=train$time
  tr_y=train$status
  te_st=test$time
  te_y=test$status
  
  set.seed(123)
  # suppressMessages(library("doParallel"))
  # registerDoParallel(detectCores())
  fit <- fit_alasso(tr_matrix, Surv(tr_st, tr_y), nfolds = 5, rule = "lambda.min", seed = c(5, 7))
  #store the model info
  lambda=fit$model$lambda
  npasses=fit$model$npasses
  devratio=fit$model$dev.ratio
  #names(fit)
  coeff=infer_variable_type(fit, tr_matrix)$name
  if(length(coeff)==0){coeff=NA}
  tr_times=train$time
  te_times=test$time
  pred_te_matrix=predict(fit,tr_matrix, Surv(tr_st, tr_y), newx = te_matrix,pred.at=te_times)
  pred_tr_matrix=predict(fit,tr_matrix, Surv(tr_st, tr_y), newx = tr_matrix,pred.at=tr_times)
  pred_te=c()
  for(i in 1:nrow(pred_te_matrix)){
    pred_te[i]=pred_te_matrix[i,][which(colnames(pred_te_matrix)==te_times[i])[1]]
  }
  pred_tr=c()
  for(i in 1:nrow(pred_tr_matrix)){
    pred_tr[i]=pred_tr_matrix[i,][which(colnames(pred_tr_matrix)==tr_times[i])[1]]
  }
  
  harrelC1 <- rcorr.cens(-pred_tr,with(train,Surv(time,status)))
  cindex_train<-harrelC1["C Index"]
  harrelC1 <- rcorr.cens(-pred_te,with(test,Surv(time,status)))
  cindex_test<-harrelC1["C Index"]
  # #calculate brier score:currently cannot be calculated if only have training data
  times=median(train$time)
  Surv.rsp <- Surv(train$time,train$status)
  Surv.rsp.new <- Surv(test$time,test$status)
  lp<-pred_tr
  lp.new=pred_te
  bs_test <- predErr(Surv.rsp,Surv.rsp.new,lp,lp.new, time = times)$error
  #maybe this is ok for the training data only (naive method)
  calculate_brier_score <- function(predicted_probs, event_status) {
    n <- length(predicted_probs)
    brier_score <- sum((predicted_probs - event_status)^2) / n
    return(brier_score)
  }
  bs_train=calculate_brier_score(lp,train$status)
  #calculate time-dependent AUC
  times=timess
  auc.uno=AUC.uno(Surv.rsp, Surv.rsp.new, lp.new, times)
  a1=auc.uno$auc[1]
  a2=auc.uno$auc[2]
  a3=auc.uno$auc[3]
  a4=auc.uno$auc[4]
  a5=auc.uno$auc[5]
  a6=auc.uno$auc[6]
  a7=auc.uno$auc[7]
  a8=auc.uno$auc[8]
  a9=auc.uno$auc[9]
  a10=auc.uno$auc[10]
  auc=auc.uno$iauc
  
  # Calculate the calibration slope:from packages rms, pec, riskregression, they all have the calibration function, however, calibration plot is dependent of the model fit, so it cannot be obtained with all models
  #calculate a naive calibration slope for the prediction
  
  #testing results on external validation datasets
  ex_matrix=data.matrix(ex_validation[,!names(ex_validation)%in% c("time","status")])
  ex_st=ex_validation$time
  ex_y=ex_validation$status
  ex_times=ex_validation$time
  pred_ex_matrix=predict(fit,ex_matrix, Surv(ex_st, ex_y), newx = ex_matrix,pred.at=ex_times)
  pred_ex=c()
  for(i in 1:nrow(pred_ex_matrix)){
    pred_ex[i]=pred_ex_matrix[i,][which(colnames(pred_ex_matrix)==ex_times[i])[1]]
  }
  harrelC1 <- rcorr.cens(-pred_ex,with(ex_validation,Surv(time,status)))
  cindex_ex<-harrelC1["C Index"]
  
  Surv.rsp <- Surv(train$time,train$status)
  Surv.rsp.new <- Surv(ex_validation$time,ex_validation$status)
  lp<-pred_tr
  lp.new=pred_ex
  times=median(train$time)
  bs_ex <- predErr(Surv.rsp,Surv.rsp.new,lp,lp.new, time = times)$error
  #calculate time-dependent AUC
  times=timess
  auc.uno=AUC.uno(Surv.rsp, Surv.rsp.new, lp.new, times)
  a1_ex=auc.uno$auc[1]
  a2_ex=auc.uno$auc[2]
  a3_ex=auc.uno$auc[3]
  a4_ex=auc.uno$auc[4]
  a5_ex=auc.uno$auc[5]
  a6_ex=auc.uno$auc[6]
  a7_ex=auc.uno$auc[7]
  a8_ex=auc.uno$auc[8]
  a9_ex=auc.uno$auc[9]
  a10_ex=auc.uno$auc[10]
  auc_ex=auc.uno$iauc
  
  return(list(coeff,cindex_test,length(coeff),bs_test,auc,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,lambda,npasses,devratio,cindex_train,bs_train,cindex_ex,bs_ex,auc_ex,a1_ex,a2_ex,a3_ex,a4_ex,a5_ex,a6_ex,a7_ex,a8_ex,a9_ex,a10_ex))
}

aenet_feature_selection=function(data,ex_validation,timess){
  test = data[[2]]
  train =data[[1]]
  tr_matrix=data.matrix(train[,!names(train)%in% c("time","status")])
  te_matrix=data.matrix(test[,!names(test)%in% c("time","status")])
  tr_st=train$time
  tr_y=train$status
  te_st=test$time
  te_y=test$status
  
  set.seed(123)
  # suppressMessages(library("doParallel"))
  # registerDoParallel(detectCores())
  fit <- fit_aenet(tr_matrix, Surv(tr_st, tr_y), nfolds = 5, rule = "lambda.min", seed = c(5, 7))
  #store the model info
  lambda=fit$model$lambda
  npasses=fit$model$npasses
  devratio=fit$model$dev.ratio
  alpha=fit$alpha
  #names(fit)
  coeff=infer_variable_type(fit, tr_matrix)$name
  if(length(coeff)==0){coeff=NA}
  tr_times=train$time
  te_times=test$time
  pred_te_matrix=predict(fit,tr_matrix, Surv(tr_st, tr_y), newx = te_matrix,pred.at=te_times)
  pred_tr_matrix=predict(fit,tr_matrix, Surv(tr_st, tr_y), newx = tr_matrix,pred.at=tr_times)
  pred_te=c()
  for(i in 1:nrow(pred_te_matrix)){
    pred_te[i]=pred_te_matrix[i,][which(colnames(pred_te_matrix)==te_times[i])[1]]
  }
  pred_tr=c()
  for(i in 1:nrow(pred_tr_matrix)){
    pred_tr[i]=pred_tr_matrix[i,][which(colnames(pred_tr_matrix)==tr_times[i])[1]]
  }
  
  harrelC1 <- rcorr.cens(-pred_tr,with(train,Surv(time,status)))
  cindex_train<-harrelC1["C Index"]
  harrelC1 <- rcorr.cens(-pred_te,with(test,Surv(time,status)))
  cindex_test<-harrelC1["C Index"]
  # #calculate brier score:currently cannot be calculated if only have training data
  times=median(train$time)
  Surv.rsp <- Surv(train$time,train$status)
  Surv.rsp.new <- Surv(test$time,test$status)
  lp<-pred_tr
  lp.new=pred_te
  bs_test <- predErr(Surv.rsp,Surv.rsp.new,lp,lp.new, time = times)$error
  #maybe this is ok for the training data only (naive method)
  calculate_brier_score <- function(predicted_probs, event_status) {
    n <- length(predicted_probs)
    brier_score <- sum((predicted_probs - event_status)^2) / n
    return(brier_score)
  }
  bs_train=calculate_brier_score(lp,train$status)
  #calculate time-dependent AUC
  times=timess
  auc.uno=AUC.uno(Surv.rsp, Surv.rsp.new, lp.new, times)
  a1=auc.uno$auc[1]
  a2=auc.uno$auc[2]
  a3=auc.uno$auc[3]
  a4=auc.uno$auc[4]
  a5=auc.uno$auc[5]
  a6=auc.uno$auc[6]
  a7=auc.uno$auc[7]
  a8=auc.uno$auc[8]
  a9=auc.uno$auc[9]
  a10=auc.uno$auc[10]
  auc=auc.uno$iauc
  
  # Calculate the calibration slope:from packages rms, pec, riskregression, they all have the calibration function, however, calibration plot is dependent of the model fit, so it cannot be obtained with all models
  #calculate a naive calibration slope for the prediction
  
  #testing results on external validation datasets
  ex_matrix=data.matrix(ex_validation[,!names(ex_validation)%in% c("time","status")])
  ex_st=ex_validation$time
  ex_y=ex_validation$status
  ex_times=ex_validation$time
  pred_ex_matrix=predict(fit,ex_matrix, Surv(ex_st, ex_y), newx = ex_matrix,pred.at=ex_times)
  pred_ex=c()
  for(i in 1:nrow(pred_ex_matrix)){
    pred_ex[i]=pred_ex_matrix[i,][which(colnames(pred_ex_matrix)==ex_times[i])[1]]
  }
  harrelC1 <- rcorr.cens(-pred_ex,with(ex_validation,Surv(time,status)))
  cindex_ex<-harrelC1["C Index"]
  
  Surv.rsp <- Surv(train$time,train$status)
  Surv.rsp.new <- Surv(ex_validation$time,ex_validation$status)
  lp<-pred_tr
  lp.new=pred_ex
  times=median(train$time)
  bs_ex <- predErr(Surv.rsp,Surv.rsp.new,lp,lp.new, time = times)$error
  #calculate time-dependent AUC
  times=timess
  auc.uno=AUC.uno(Surv.rsp, Surv.rsp.new, lp.new, times)
  a1_ex=auc.uno$auc[1]
  a2_ex=auc.uno$auc[2]
  a3_ex=auc.uno$auc[3]
  a4_ex=auc.uno$auc[4]
  a5_ex=auc.uno$auc[5]
  a6_ex=auc.uno$auc[6]
  a7_ex=auc.uno$auc[7]
  a8_ex=auc.uno$auc[8]
  a9_ex=auc.uno$auc[9]
  a10_ex=auc.uno$auc[10]
  auc_ex=auc.uno$iauc
  
  return(list(coeff,cindex_test,length(coeff),bs_test,auc,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,lambda,npasses,devratio,cindex_train,bs_train,cindex_ex,bs_ex,auc_ex,a1_ex,a2_ex,a3_ex,a4_ex,a5_ex,a6_ex,a7_ex,a8_ex,a9_ex,a10_ex,alpha))
}

mcp_feature_selection=function(data,ex_validation,timess){
  test = data[[2]]
  train =data[[1]]
  tr_matrix=data.matrix(train[,!names(train)%in% c("time","status")])
  te_matrix=data.matrix(test[,!names(test)%in% c("time","status")])
  tr_st=train$time
  tr_y=train$status
  te_st=test$time
  te_y=test$status
  
  set.seed(123)
  # suppressMessages(library("doParallel"))
  # registerDoParallel(detectCores())
  fit <- fit_mcp(tr_matrix, Surv(tr_st, tr_y), nfolds = 5)
  #store the model info
  lambda=fit$model$lambda
  gamma=fit$model$gamma
  iter=fit$model$iter
  loss=fit$model$loss
  #names(fit)
  coeff=infer_variable_type(fit, tr_matrix)$name
  if(length(coeff)==0){coeff=NA}
  tr_times=train$time
  te_times=test$time
  pred_te_matrix=predict(fit,tr_matrix, Surv(tr_st, tr_y), newx = te_matrix,pred.at=te_times)
  pred_tr_matrix=predict(fit,tr_matrix, Surv(tr_st, tr_y), newx = tr_matrix,pred.at=tr_times)
  pred_te=c()
  for(i in 1:nrow(pred_te_matrix)){
    pred_te[i]=pred_te_matrix[i,][which(colnames(pred_te_matrix)==te_times[i])[1]]
  }
  pred_tr=c()
  for(i in 1:nrow(pred_tr_matrix)){
    pred_tr[i]=pred_tr_matrix[i,][which(colnames(pred_tr_matrix)==tr_times[i])[1]]
  }
  
  harrelC1 <- rcorr.cens(-pred_tr,with(train,Surv(time,status)))
  cindex_train<-harrelC1["C Index"]
  harrelC1 <- rcorr.cens(-pred_te,with(test,Surv(time,status)))
  cindex_test<-harrelC1["C Index"]
  # #calculate brier score:currently cannot be calculated if only have training data
  times=median(train$time)
  Surv.rsp <- Surv(train$time,train$status)
  Surv.rsp.new <- Surv(test$time,test$status)
  lp<-pred_tr
  lp.new=pred_te
  bs_test <- predErr(Surv.rsp,Surv.rsp.new,lp,lp.new, time = times)$error
  #maybe this is ok for the training data only (naive method)
  calculate_brier_score <- function(predicted_probs, event_status) {
    n <- length(predicted_probs)
    brier_score <- sum((predicted_probs - event_status)^2) / n
    return(brier_score)
  }
  bs_train=calculate_brier_score(lp,train$status)
  #calculate time-dependent AUC
  times=timess
  auc.uno=AUC.uno(Surv.rsp, Surv.rsp.new, lp.new, times)
  a1=auc.uno$auc[1]
  a2=auc.uno$auc[2]
  a3=auc.uno$auc[3]
  a4=auc.uno$auc[4]
  a5=auc.uno$auc[5]
  a6=auc.uno$auc[6]
  a7=auc.uno$auc[7]
  a8=auc.uno$auc[8]
  a9=auc.uno$auc[9]
  a10=auc.uno$auc[10]
  auc=auc.uno$iauc
  
  # Calculate the calibration slope:from packages rms, pec, riskregression, they all have the calibration function, however, calibration plot is dependent of the model fit, so it cannot be obtained with all models
  #calculate a naive calibration slope for the prediction
  
  #testing results on external validation datasets
  ex_matrix=data.matrix(ex_validation[,!names(ex_validation)%in% c("time","status")])
  ex_st=ex_validation$time
  ex_y=ex_validation$status
  ex_times=ex_validation$time
  pred_ex_matrix=predict(fit,ex_matrix, Surv(ex_st, ex_y), newx = ex_matrix,pred.at=ex_times)
  pred_ex=c()
  for(i in 1:nrow(pred_ex_matrix)){
    pred_ex[i]=pred_ex_matrix[i,][which(colnames(pred_ex_matrix)==ex_times[i])[1]]
  }
  harrelC1 <- rcorr.cens(-pred_ex,with(ex_validation,Surv(time,status)))
  cindex_ex<-harrelC1["C Index"]
  
  Surv.rsp <- Surv(train$time,train$status)
  Surv.rsp.new <- Surv(ex_validation$time,ex_validation$status)
  lp<-pred_tr
  lp.new=pred_ex
  times=median(train$time)
  bs_ex <- predErr(Surv.rsp,Surv.rsp.new,lp,lp.new, time = times)$error
  #calculate time-dependent AUC
  times=timess
  auc.uno=AUC.uno(Surv.rsp, Surv.rsp.new, lp.new, times)
  a1_ex=auc.uno$auc[1]
  a2_ex=auc.uno$auc[2]
  a3_ex=auc.uno$auc[3]
  a4_ex=auc.uno$auc[4]
  a5_ex=auc.uno$auc[5]
  a6_ex=auc.uno$auc[6]
  a7_ex=auc.uno$auc[7]
  a8_ex=auc.uno$auc[8]
  a9_ex=auc.uno$auc[9]
  a10_ex=auc.uno$auc[10]
  auc_ex=auc.uno$iauc
  
  return(list(coeff,cindex_test,length(coeff),bs_test,auc,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,lambda,gamma,iter,cindex_train,bs_train,cindex_ex,bs_ex,auc_ex,a1_ex,a2_ex,a3_ex,a4_ex,a5_ex,a6_ex,a7_ex,a8_ex,a9_ex,a10_ex,loss))
}

mnet_feature_selection=function(data,ex_validation,timess){
  test = data[[2]]
  train =data[[1]]
  tr_matrix=data.matrix(train[,!names(train)%in% c("time","status")])
  te_matrix=data.matrix(test[,!names(test)%in% c("time","status")])
  tr_st=train$time
  tr_y=train$status
  te_st=test$time
  te_y=test$status
  
  set.seed(123)
  # suppressMessages(library("doParallel"))
  # registerDoParallel(detectCores())
  fit <- fit_mnet(tr_matrix, Surv(tr_st, tr_y), nfolds = 5)
  #store the model info
  lambda=fit$model$lambda
  gamma=fit$model$gamma
  iter=fit$model$iter
  loss=fit$model$loss
  alpha=fit$model$alpha
  #names(fit)
  coeff=infer_variable_type(fit, tr_matrix)$name
  if(length(coeff)==0){coeff=NA}
  tr_times=train$time
  te_times=test$time
  pred_te_matrix=predict(fit,tr_matrix, Surv(tr_st, tr_y), newx = te_matrix,pred.at=te_times)
  pred_tr_matrix=predict(fit,tr_matrix, Surv(tr_st, tr_y), newx = tr_matrix,pred.at=tr_times)
  pred_te=c()
  for(i in 1:nrow(pred_te_matrix)){
    pred_te[i]=pred_te_matrix[i,][which(colnames(pred_te_matrix)==te_times[i])[1]]
  }
  pred_tr=c()
  for(i in 1:nrow(pred_tr_matrix)){
    pred_tr[i]=pred_tr_matrix[i,][which(colnames(pred_tr_matrix)==tr_times[i])[1]]
  }
  
  harrelC1 <- rcorr.cens(-pred_tr,with(train,Surv(time,status)))
  cindex_train<-harrelC1["C Index"]
  harrelC1 <- rcorr.cens(-pred_te,with(test,Surv(time,status)))
  cindex_test<-harrelC1["C Index"]
  # #calculate brier score:currently cannot be calculated if only have training data
  times=median(train$time)
  Surv.rsp <- Surv(train$time,train$status)
  Surv.rsp.new <- Surv(test$time,test$status)
  lp<-pred_tr
  lp.new=pred_te
  bs_test <- predErr(Surv.rsp,Surv.rsp.new,lp,lp.new, time = times)$error
  #maybe this is ok for the training data only (naive method)
  calculate_brier_score <- function(predicted_probs, event_status) {
    n <- length(predicted_probs)
    brier_score <- sum((predicted_probs - event_status)^2) / n
    return(brier_score)
  }
  bs_train=calculate_brier_score(lp,train$status)
  #calculate time-dependent AUC
  times=timess
  auc.uno=AUC.uno(Surv.rsp, Surv.rsp.new, lp.new, times)
  a1=auc.uno$auc[1]
  a2=auc.uno$auc[2]
  a3=auc.uno$auc[3]
  a4=auc.uno$auc[4]
  a5=auc.uno$auc[5]
  a6=auc.uno$auc[6]
  a7=auc.uno$auc[7]
  a8=auc.uno$auc[8]
  a9=auc.uno$auc[9]
  a10=auc.uno$auc[10]
  auc=auc.uno$iauc
  
  # Calculate the calibration slope:from packages rms, pec, riskregression, they all have the calibration function, however, calibration plot is dependent of the model fit, so it cannot be obtained with all models
  #calculate a naive calibration slope for the prediction
  
  #testing results on external validation datasets
  ex_matrix=data.matrix(ex_validation[,!names(ex_validation)%in% c("time","status")])
  ex_st=ex_validation$time
  ex_y=ex_validation$status
  ex_times=ex_validation$time
  pred_ex_matrix=predict(fit,ex_matrix, Surv(ex_st, ex_y), newx = ex_matrix,pred.at=ex_times)
  pred_ex=c()
  for(i in 1:nrow(pred_ex_matrix)){
    pred_ex[i]=pred_ex_matrix[i,][which(colnames(pred_ex_matrix)==ex_times[i])[1]]
  }
  harrelC1 <- rcorr.cens(-pred_ex,with(ex_validation,Surv(time,status)))
  cindex_ex<-harrelC1["C Index"]
  
  Surv.rsp <- Surv(train$time,train$status)
  Surv.rsp.new <- Surv(ex_validation$time,ex_validation$status)
  lp<-pred_tr
  lp.new=pred_ex
  times=median(train$time)
  bs_ex <- predErr(Surv.rsp,Surv.rsp.new,lp,lp.new, time = times)$error
  #calculate time-dependent AUC
  times=timess
  auc.uno=AUC.uno(Surv.rsp, Surv.rsp.new, lp.new, times)
  a1_ex=auc.uno$auc[1]
  a2_ex=auc.uno$auc[2]
  a3_ex=auc.uno$auc[3]
  a4_ex=auc.uno$auc[4]
  a5_ex=auc.uno$auc[5]
  a6_ex=auc.uno$auc[6]
  a7_ex=auc.uno$auc[7]
  a8_ex=auc.uno$auc[8]
  a9_ex=auc.uno$auc[9]
  a10_ex=auc.uno$auc[10]
  auc_ex=auc.uno$iauc
  
  return(list(coeff,cindex_test,length(coeff),bs_test,auc,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,lambda,gamma,iter,cindex_train,bs_train,cindex_ex,bs_ex,auc_ex,a1_ex,a2_ex,a3_ex,a4_ex,a5_ex,a6_ex,a7_ex,a8_ex,a9_ex,a10_ex,loss,alpha))
}

scad_feature_selection=function(data,ex_validation,timess){
  test = data[[2]]
  train =data[[1]]
  tr_matrix=data.matrix(train[,!names(train)%in% c("time","status")])
  te_matrix=data.matrix(test[,!names(test)%in% c("time","status")])
  tr_st=train$time
  tr_y=train$status
  te_st=test$time
  te_y=test$status
  
  set.seed(123)
  # suppressMessages(library("doParallel"))
  # registerDoParallel(detectCores())
  fit <- fit_scad(tr_matrix, Surv(tr_st, tr_y), nfolds = 5)
  #store the model info
  lambda=fit$model$lambda
  gamma=fit$model$gamma
  iter=fit$model$iter
  loss=fit$model$loss
  #names(fit)
  coeff=infer_variable_type(fit, tr_matrix)$name
  if(length(coeff)==0){coeff=NA}
  tr_times=train$time
  te_times=test$time
  pred_te_matrix=predict(fit,tr_matrix, Surv(tr_st, tr_y), newx = te_matrix,pred.at=te_times)
  pred_tr_matrix=predict(fit,tr_matrix, Surv(tr_st, tr_y), newx = tr_matrix,pred.at=tr_times)
  pred_te=c()
  for(i in 1:nrow(pred_te_matrix)){
    pred_te[i]=pred_te_matrix[i,][which(colnames(pred_te_matrix)==te_times[i])[1]]
  }
  pred_tr=c()
  for(i in 1:nrow(pred_tr_matrix)){
    pred_tr[i]=pred_tr_matrix[i,][which(colnames(pred_tr_matrix)==tr_times[i])[1]]
  }
  
  harrelC1 <- rcorr.cens(-pred_tr,with(train,Surv(time,status)))
  cindex_train<-harrelC1["C Index"]
  harrelC1 <- rcorr.cens(-pred_te,with(test,Surv(time,status)))
  cindex_test<-harrelC1["C Index"]
  # #calculate brier score:currently cannot be calculated if only have training data
  times=median(train$time)
  Surv.rsp <- Surv(train$time,train$status)
  Surv.rsp.new <- Surv(test$time,test$status)
  lp<-pred_tr
  lp.new=pred_te
  bs_test <- predErr(Surv.rsp,Surv.rsp.new,lp,lp.new, time = times)$error
  #maybe this is ok for the training data only (naive method)
  calculate_brier_score <- function(predicted_probs, event_status) {
    n <- length(predicted_probs)
    brier_score <- sum((predicted_probs - event_status)^2) / n
    return(brier_score)
  }
  bs_train=calculate_brier_score(lp,train$status)
  #calculate time-dependent AUC
  times=timess
  auc.uno=AUC.uno(Surv.rsp, Surv.rsp.new, lp.new, times)
  a1=auc.uno$auc[1]
  a2=auc.uno$auc[2]
  a3=auc.uno$auc[3]
  a4=auc.uno$auc[4]
  a5=auc.uno$auc[5]
  a6=auc.uno$auc[6]
  a7=auc.uno$auc[7]
  a8=auc.uno$auc[8]
  a9=auc.uno$auc[9]
  a10=auc.uno$auc[10]
  auc=auc.uno$iauc
  
  # Calculate the calibration slope:from packages rms, pec, riskregression, they all have the calibration function, however, calibration plot is dependent of the model fit, so it cannot be obtained with all models
  #calculate a naive calibration slope for the prediction
  
  #testing results on external validation datasets
  ex_matrix=data.matrix(ex_validation[,!names(ex_validation)%in% c("time","status")])
  ex_st=ex_validation$time
  ex_y=ex_validation$status
  ex_times=ex_validation$time
  pred_ex_matrix=predict(fit,ex_matrix, Surv(ex_st, ex_y), newx = ex_matrix,pred.at=ex_times)
  pred_ex=c()
  for(i in 1:nrow(pred_ex_matrix)){
    pred_ex[i]=pred_ex_matrix[i,][which(colnames(pred_ex_matrix)==ex_times[i])[1]]
  }
  harrelC1 <- rcorr.cens(-pred_ex,with(ex_validation,Surv(time,status)))
  cindex_ex<-harrelC1["C Index"]
  
  Surv.rsp <- Surv(train$time,train$status)
  Surv.rsp.new <- Surv(ex_validation$time,ex_validation$status)
  lp<-pred_tr
  lp.new=pred_ex
  times=median(train$time)
  bs_ex <- predErr(Surv.rsp,Surv.rsp.new,lp,lp.new, time = times)$error
  #calculate time-dependent AUC
  times=timess
  auc.uno=AUC.uno(Surv.rsp, Surv.rsp.new, lp.new, times)
  a1_ex=auc.uno$auc[1]
  a2_ex=auc.uno$auc[2]
  a3_ex=auc.uno$auc[3]
  a4_ex=auc.uno$auc[4]
  a5_ex=auc.uno$auc[5]
  a6_ex=auc.uno$auc[6]
  a7_ex=auc.uno$auc[7]
  a8_ex=auc.uno$auc[8]
  a9_ex=auc.uno$auc[9]
  a10_ex=auc.uno$auc[10]
  auc_ex=auc.uno$iauc
  
  return(list(coeff,cindex_test,length(coeff),bs_test,auc,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,lambda,gamma,iter,cindex_train,bs_train,cindex_ex,bs_ex,auc_ex,a1_ex,a2_ex,a3_ex,a4_ex,a5_ex,a6_ex,a7_ex,a8_ex,a9_ex,a10_ex,loss))
}

snet_feature_selection=function(data,ex_validation,timess){
  test = data[[2]]
  train =data[[1]]
  tr_matrix=data.matrix(train[,!names(train)%in% c("time","status")])
  te_matrix=data.matrix(test[,!names(test)%in% c("time","status")])
  tr_st=train$time
  tr_y=train$status
  te_st=test$time
  te_y=test$status
  
  set.seed(123)
  # suppressMessages(library("doParallel"))
  # registerDoParallel(detectCores())
  fit <- fit_snet(tr_matrix, Surv(tr_st, tr_y), nfolds = 5)
  #store the model info
  lambda=fit$model$lambda
  gamma=fit$model$gamma
  iter=fit$model$iter
  loss=fit$model$loss
  alpha=fit$model$alpha
  #names(fit)
  coeff=infer_variable_type(fit, tr_matrix)$name
  if(length(coeff)==0){coeff=NA}
  tr_times=train$time
  te_times=test$time
  pred_te_matrix=predict(fit,tr_matrix, Surv(tr_st, tr_y), newx = te_matrix,pred.at=te_times)
  pred_tr_matrix=predict(fit,tr_matrix, Surv(tr_st, tr_y), newx = tr_matrix,pred.at=tr_times)
  pred_te=c()
  for(i in 1:nrow(pred_te_matrix)){
    pred_te[i]=pred_te_matrix[i,][which(colnames(pred_te_matrix)==te_times[i])[1]]
  }
  pred_tr=c()
  for(i in 1:nrow(pred_tr_matrix)){
    pred_tr[i]=pred_tr_matrix[i,][which(colnames(pred_tr_matrix)==tr_times[i])[1]]
  }
  
  harrelC1 <- rcorr.cens(-pred_tr,with(train,Surv(time,status)))
  cindex_train<-harrelC1["C Index"]
  harrelC1 <- rcorr.cens(-pred_te,with(test,Surv(time,status)))
  cindex_test<-harrelC1["C Index"]
  # #calculate brier score:currently cannot be calculated if only have training data
  times=median(train$time)
  Surv.rsp <- Surv(train$time,train$status)
  Surv.rsp.new <- Surv(test$time,test$status)
  lp<-pred_tr
  lp.new=pred_te
  bs_test <- predErr(Surv.rsp,Surv.rsp.new,lp,lp.new, time = times)$error
  #maybe this is ok for the training data only (naive method)
  calculate_brier_score <- function(predicted_probs, event_status) {
    n <- length(predicted_probs)
    brier_score <- sum((predicted_probs - event_status)^2) / n
    return(brier_score)
  }
  bs_train=calculate_brier_score(lp,train$status)
  #calculate time-dependent AUC
  times=timess
  auc.uno=AUC.uno(Surv.rsp, Surv.rsp.new, lp.new, times)
  a1=auc.uno$auc[1]
  a2=auc.uno$auc[2]
  a3=auc.uno$auc[3]
  a4=auc.uno$auc[4]
  a5=auc.uno$auc[5]
  a6=auc.uno$auc[6]
  a7=auc.uno$auc[7]
  a8=auc.uno$auc[8]
  a9=auc.uno$auc[9]
  a10=auc.uno$auc[10]
  auc=auc.uno$iauc
  
  # Calculate the calibration slope:from packages rms, pec, riskregression, they all have the calibration function, however, calibration plot is dependent of the model fit, so it cannot be obtained with all models
  #calculate a naive calibration slope for the prediction
  
  #testing results on external validation datasets
  ex_matrix=data.matrix(ex_validation[,!names(ex_validation)%in% c("time","status")])
  ex_st=ex_validation$time
  ex_y=ex_validation$status
  ex_times=ex_validation$time
  pred_ex_matrix=predict(fit,ex_matrix, Surv(ex_st, ex_y), newx = ex_matrix,pred.at=ex_times)
  pred_ex=c()
  for(i in 1:nrow(pred_ex_matrix)){
    pred_ex[i]=pred_ex_matrix[i,][which(colnames(pred_ex_matrix)==ex_times[i])[1]]
  }
  harrelC1 <- rcorr.cens(-pred_ex,with(ex_validation,Surv(time,status)))
  cindex_ex<-harrelC1["C Index"]
  
  Surv.rsp <- Surv(train$time,train$status)
  Surv.rsp.new <- Surv(ex_validation$time,ex_validation$status)
  lp<-pred_tr
  lp.new=pred_ex
  times=median(train$time)
  bs_ex <- predErr(Surv.rsp,Surv.rsp.new,lp,lp.new, time = times)$error
  #calculate time-dependent AUC
  times=timess
  auc.uno=AUC.uno(Surv.rsp, Surv.rsp.new, lp.new, times)
  a1_ex=auc.uno$auc[1]
  a2_ex=auc.uno$auc[2]
  a3_ex=auc.uno$auc[3]
  a4_ex=auc.uno$auc[4]
  a5_ex=auc.uno$auc[5]
  a6_ex=auc.uno$auc[6]
  a7_ex=auc.uno$auc[7]
  a8_ex=auc.uno$auc[8]
  a9_ex=auc.uno$auc[9]
  a10_ex=auc.uno$auc[10]
  auc_ex=auc.uno$iauc
  
  return(list(coeff,cindex_test,length(coeff),bs_test,auc,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,lambda,gamma,iter,cindex_train,bs_train,cindex_ex,bs_ex,auc_ex,a1_ex,a2_ex,a3_ex,a4_ex,a5_ex,a6_ex,a7_ex,a8_ex,a9_ex,a10_ex,loss,alpha))
}