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

#############################################################
library(glmnet)
library(PRROC)
# predictor : logistic classifier from training set
wts <- rep(1,length(ytrain))
wts[which(ytrain==1)] <- length(ytrain)/(2*length(which(ytrain==1))) #14.0231214
wts[which(ytrain!=1)] <- length(ytrain)/(2*length(which(ytrain!=1))) #0.5184869

# lasso
ftr <- ftr_sel(xtrain,100000)
Xtrain_sub <- xtrain[,ftr]
set.seed(123)
fit <- cv.glmnet(Xtrain_sub, ytrain, alpha=1, family="binomial", nfolds=3, weights=wts)


# enet
set.seed(123)
fit2 <- cv.glmnet(xtrain, ytrain, alpha=0.5, family="binomial", nfolds=3, weights=wts)


# ridge
ftr <- ftr_sel(xtrain,300000)
Xtrain_sub <- xtrain[,ftr]
set.seed(123)
fit3 <- cv.glmnet(Xtrain_sub, ytrain, alpha=0, family="binomial", nfolds=3, weights=wts)


######################################################################

logistic_predictor <- predict(fit,xtest[,ftr],type="link") #lasso
# logistic_predictor <- predict(fit2,xtest[,ftr],type="link") #enet
# logistic_predictor <- predict(fit3,xtest[,ftr],type="link") #ridge
logistic_predictor_data <- cbind.data.frame(time=age_event_test,status=ytest,logistic_predictor,age=agetest,sex=sextest)


#########################################################################
cox_model<- summary(coxph(Surv(time,status) ~ scale(logistic_predictor)+age+sex, logistic_predictor_data))
pseudo_r2 <- cox_model$rsq[[1]]
HR <- cox_model$conf.int[,-2]%>%cbind(pval=cox_model$coefficients[,5])
coef <- coef(fit,s="lambda.min")
active_idx <- which(coef!=0)
n_features <- length(active_idx)-1
cox_output <- cbind.data.frame(predictor="Logistic classifier from the training set",
                               COV=300000,
                               n_features=n_features,
                               pseudo_r2=pseudo_r2,
                               X=rownames(HR),
                               HR)
