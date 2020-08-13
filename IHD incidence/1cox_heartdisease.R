library(doMC)
library(tidyverse)
library(glmnet)
library(survival)


# create a function to calculate COV for all CpGs and store the CpG ids for the thresholds specified in the n_ftr argument
ftr_sel <- function(x,n_ftr) {
  sd=apply(x,2,sd)
  mean=apply(x,2,mean)
  cov = abs(sd/mean)
  sorted_cov <- sort(cov,index.return=TRUE,decreasing=TRUE)
  CpG_id <- colnames(x)[sorted_cov$ix]
  CpG_id[1:n_ftr]
}

###########################################################################

pip_cox <- function(xtrain,ytrain,age_event_train,xtest,ytest,age_event_test,penalty,outcome,nfolds_cv,COV) {
  pseudo_r2 <- numeric(length(COV))
  
  # set class weights
  wts <- rep(1,length(ytrain))
  wts[which(ytrain==1)] <- length(ytrain)/(2*length(which(ytrain==1))) #14.0231214
  wts[which(ytrain!=1)] <- length(ytrain)/(2*length(which(ytrain!=1))) #0.5184869
  
  for (k in 1:length(COV)) {
    ftr <- ftr_sel(xtrain,COV[k])
    # Construct X training subset
    Xtrain_sub <- xtrain[,ftr]
    
    registerDoMC(cores=2)
    set.seed(123)
    fit <- cv.glmnet(Xtrain_sub,Surv(age_event_train,ytrain),alpha=penalty,family=outcome,nfolds=nfolds_cv,weights=wts,parallel=T)
    # extract coefficients
    coef <- coef(fit, s = "lambda.min")
    active_coef <- as.data.frame(coef[which(coef!=0),])
    write_csv(active_coef,path=paste0("/Cluster_Filespace/Marioni_Group/Yufei/output/cox/coef/",penalty,"coef",COV[k],"_cox.csv"))
    n_features <- nrow(active_coef)
    
    # metrics on the training set
    fit_predictor <- predict(fit,Xtrain_sub,type="link",s="lambda.min")[,1]
    fit_predictor_data <- cbind.data.frame(time=age_event_train,status=ytrain,fit_predictor=fit_predictor,age=agetrain,sex=sextrain)
    fit_model <- summary(coxph(Surv(time,status) ~ scale(fit_predictor)+age+sex,fit_predictor_data))
    pseudo_r2[k] <- fit_model$rsq[[1]]
    HR <- fit_model$conf.int[,-2]%>%cbind(pval=fit_model$coefficients[,5])
    fit_output <- cbind.data.frame(penalty=penalty,
                                   dataset="training set",
                                   COV=COV[k],
                                   n_features=n_features,
                                   pseudo_r2=pseudo_r2[k],
                                   X=rownames(HR),
                                   HR)
    write_csv(fit_output,path=paste0("/Cluster_Filespace/Marioni_Group/Yufei/output/cox/coxpredictortrain/",penalty,"coxpredictor_withagesex",COV[k],"_cox.csv"))
    
    # metrics on the test set
    cox_predictor <- predict(fit,xtest[,ftr],type="link",s="lambda.min")[,1]
    cox_predictor_data <- cbind.data.frame(time=age_event_test,status=ytest,cox_predictor=cox_predictor,age=agetest,sex=sextest)
    cox_model2<- summary(coxph(Surv(time,status) ~ scale(cox_predictor)+age+sex, cox_predictor_data))
    pseudo_r22 <- cox_model2$rsq[[1]]
    HR2 <- cox_model2$conf.int[,-2]%>%cbind(pval = cox_model2$coefficients[,5])
    cox_output2 <- cbind.data.frame(penalty=penalty,
                                    predictor="cox predictor from the training set",
                                    COV=COV[k],
                                    n_features=n_features,
                                    pseudo_r22=pseudo_r22,
                                    X=rownames(HR2),
                                    HR2)
    write_csv(cox_output2, path=paste0("/Cluster_Filespace/Marioni_Group/Yufei/output/cox/coxpredictortest/",penalty,"coxpredictor_withagesex",COV[k],"_cox.csv"))
  }
  # find the best model
  m <- which(pseudo_r2==max(pseudo_r2))[1]
  cat("When penalty = ",penalty," the best model: COV = ", COV[m])
}

COV <- c(10000,20000,30000,40000,50000,100000,200000,300000,ncol(xtrain))
COV <- ncol(xtrain)
system.time(
  pip_cox(xtrain,ytrain,age_event_train,xtest,ytest,age_event_test,penalty=1,outcome="cox",nfolds_cv=3,COV)
)
system.time(
  pip_cox(xtrain,ytrain,age_event_train,xtest,ytest,age_event_test,penalty=0,outcome="cox",nfolds_cv=3,COV)
)
system.time(
  pip_cox(xtrain,ytrain,age_event_train,xtest,ytest,age_event_test,penalty=0.5,outcome="cox",nfolds_cv=3,COV)
)
