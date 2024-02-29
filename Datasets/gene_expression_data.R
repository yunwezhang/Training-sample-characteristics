# this file is used to obtain curated gene expression datasets with clinical information gender based on the paper: csurvival
# https://tau.cmmt.ubc.ca/cSurvival/

#ACC: Adrenocortical Carcinoma, female:male: 1.4:1 (more in female: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7036530/ ref22)

setwd("/dskh/nobackup/yunwei/stability_selection_codes")
acc_clinical=read.csv("/dskh/nobackup/yunwei/stability_selection_codes/TCGA-ACC_clinical.csv")
acc_clinical=acc_clinical[,c("patient_id","gender","OS","OS.time")]
dim(acc_clinical)

acc_expression=read.csv("/dskh/nobackup/yunwei/stability_selection_codes/TCGA-ACC_expression.csv",header = TRUE)
dim(acc_expression)

acc_clinical2=acc_clinical[match(acc_expression$patient_id,acc_clinical$patient_id),]
acc_expression2=cbind.data.frame(acc_expression,acc_clinical2)
dim(acc_expression2)
acc_expression3=acc_expression2[,-which(colnames(acc_expression2)%in%c("patient_id"))]
dim(acc_expression3)
#get rid of 0 entries all the way
not_all_na <- function(x) any(!x==0)
check=apply(acc_expression3,2,not_all_na)
check1=as.data.frame(check)
any(check1$check==FALSE)
sum(check1$check==FALSE)
names=rownames(check1)[which(check1$check==FALSE)]
acc_expression4=acc_expression3[,which(!colnames(acc_expression3)%in%names)]
dim(acc_expression4)

#log2 transform
acc_expression5=acc_expression4
acc_expression5[,1:(dim(acc_expression4)[2]-3)]=log2(acc_expression4[,1:(dim(acc_expression4)[2]-3)]+1)
#quantile normalisation
dim(acc_expression5)
trans_matrix=t(acc_expression5[,1:(dim(acc_expression5)[2]-3)])
trans_matrix[1:6,1:6]
library(preprocessCore)
data_norm <- normalize.quantiles(trans_matrix, copy = TRUE)
data_norm[1:6,1:6]
boxplot(data_norm)
norm_matrix=as.data.frame(t(data_norm))
colnames(norm_matrix)=rownames(trans_matrix)
acc_expression6=cbind.data.frame(norm_matrix,acc_expression5[,(dim(acc_expression5)[2]-2):dim(acc_expression5)[2]])
acc_expression6[1:6,1:6]
write.csv(acc_expression6,file = "ACC.csv")


################

setwd("/dskh/nobackup/yunwei/stability_selection_codes")
acc_clinical=read.csv("/dskh/nobackup/yunwei/stability_selection_codes/TCGA-BRCA_clinical.csv")
acc_clinical=acc_clinical[,c("patient_id","gender","OS","OS.time")]
dim(acc_clinical)

acc_expression=read.csv("/dskh/nobackup/yunwei/stability_selection_codes/TCGA-BRCA_expression.csv",header = TRUE)
dim(acc_expression)

acc_clinical2=acc_clinical[match(acc_expression$patient_id,acc_clinical$patient_id),]
acc_expression2=cbind.data.frame(acc_expression,acc_clinical2)
dim(acc_expression2)
acc_expression3=acc_expression2[,-which(colnames(acc_expression2)%in%c("patient_id"))]
dim(acc_expression3)
#get rid of 0 entries all the way
not_all_na <- function(x) any(!x==0)
check=apply(acc_expression3,2,not_all_na)
check1=as.data.frame(check)
any(check1$check==FALSE)
sum(check1$check==FALSE)
names=rownames(check1)[which(check1$check==FALSE)]
acc_expression4=acc_expression3[,which(!colnames(acc_expression3)%in%names)]
dim(acc_expression4)

#log2 transform
acc_expression5=acc_expression4
acc_expression5[,1:(dim(acc_expression4)[2]-3)]=log2(acc_expression4[,1:(dim(acc_expression4)[2]-3)]+1)
#quantile normalisation
dim(acc_expression5)
trans_matrix=t(acc_expression5[,1:(dim(acc_expression5)[2]-3)])
trans_matrix[1:6,1:6]
library(preprocessCore)
data_norm <- normalize.quantiles(trans_matrix, copy = TRUE)
data_norm[1:6,1:6]
boxplot(data_norm)
norm_matrix=as.data.frame(t(data_norm))
colnames(norm_matrix)=rownames(trans_matrix)
acc_expression6=cbind.data.frame(norm_matrix,acc_expression5[,(dim(acc_expression5)[2]-2):dim(acc_expression5)[2]])
acc_expression6[1:6,1:6]
write.csv(acc_expression6,file = "BRCA.csv")


################

setwd("/dskh/nobackup/yunwei/stability_selection_codes")
acc_clinical=read.csv("/dskh/nobackup/yunwei/stability_selection_codes/TCGA-BLCA_clinical.csv")
acc_clinical=acc_clinical[,c("patient_id","gender","OS","OS.time")]
dim(acc_clinical)

acc_expression=read.csv("/dskh/nobackup/yunwei/stability_selection_codes/TCGA-BLCA_expression.csv",header = TRUE)
dim(acc_expression)

acc_clinical2=acc_clinical[match(acc_expression$patient_id,acc_clinical$patient_id),]
acc_expression2=cbind.data.frame(acc_expression,acc_clinical2)
dim(acc_expression2)
acc_expression3=acc_expression2[,-which(colnames(acc_expression2)%in%c("patient_id"))]
dim(acc_expression3)
#get rid of 0 entries all the way
not_all_na <- function(x) any(!x==0)
check=apply(acc_expression3,2,not_all_na)
check1=as.data.frame(check)
any(check1$check==FALSE)
sum(check1$check==FALSE)
names=rownames(check1)[which(check1$check==FALSE)]
acc_expression4=acc_expression3[,which(!colnames(acc_expression3)%in%names)]
dim(acc_expression4)

#log2 transform
acc_expression5=acc_expression4
acc_expression5[,1:(dim(acc_expression4)[2]-3)]=log2(acc_expression4[,1:(dim(acc_expression4)[2]-3)]+1)
#quantile normalisation
dim(acc_expression5)
trans_matrix=t(acc_expression5[,1:(dim(acc_expression5)[2]-3)])
trans_matrix[1:6,1:6]
library(preprocessCore)
data_norm <- normalize.quantiles(trans_matrix, copy = TRUE)
data_norm[1:6,1:6]
boxplot(data_norm)
norm_matrix=as.data.frame(t(data_norm))
colnames(norm_matrix)=rownames(trans_matrix)
acc_expression6=cbind.data.frame(norm_matrix,acc_expression5[,(dim(acc_expression5)[2]-2):dim(acc_expression5)[2]])
acc_expression6[1:6,1:6]
write.csv(acc_expression6,file = "BLCA.csv")

