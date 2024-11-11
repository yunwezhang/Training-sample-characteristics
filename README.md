# Training-sample-characteristics

This repository provides codes for examining how training sample characteristics contribute to survival model performances as an empirical study. We provide codes to generate figures for our manuscript along with an example code file to get results for a real-life dataset. As codes will take more than an hour to run on a single-core laptop/pc. We recommend users to use our saved .rds and .rdata files. All real datasets used in our study are publicly available and please email the corresponding author to obtain. 

## This repo 
We have our high-resolution figures under the Figures folder. Codes for results under the Codes folder. Example for one real dataset under the Example folder.

## Package info
We require the following packages to work to run our codes.
<details>
<summary>Click to expand</summary>
  
```r
library(survAUC)
library(reshape2)
library(ggplot2)
library(readxl)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(MASS)
library(survival)
library(DT)
library(EnvStats)
library(pbmcapply)
library(pec)
library(survivalROC)
library(survival)
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
library(hdnom)
library(Hmisc)
library(ggrepel)

```

