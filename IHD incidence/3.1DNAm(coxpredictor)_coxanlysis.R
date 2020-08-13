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
library(doMC)
ftr <- ftr_sel(xtrain,ncol(xtrain))

Xtrain_sub <- xtrain[,ftr]

# set class weights
wts <- rep(1,length(ytrain))
wts[which(ytrain==1)] <- length(ytrain)/(2*length(which(ytrain==1))) #14.0231214
wts[which(ytrain!=1)] <- length(ytrain)/(2*length(which(ytrain!=1))) #0.5184869

penalty <- 1 # or 0.5 or 0 (lasso/enet/ridge)
registerDoMC(cores=2)
set.seed(123)
fit <- cv.glmnet(Xtrain_sub,Surv(age_event_train,ytrain),alpha=penalty,family="cox",nfolds=3,weights=wts,parallel=T)
# extract coefficients
coef <- coef(fit, s = "lambda.min")
active_idx <- which(coef!=0)
n_features <- length(active_idx)


cox_predictor <- predict(fit,xtest[,ftr],type="link",s="lambda.min")[,1]
cox_predictor_data <- cbind.data.frame(time=age_event_test,status=ytest,cox_predictor=cox_predictor,age=agetest,sex=sextest)

#########################################################################

cox_model<- summary(coxph(Surv(time,status) ~ scale(cox_predictor)+age+sex, cox_predictor_data))
pseudo_r2 <- cox_model$rsq[[1]]
HR <- cox_model$conf.int[,-2]%>%cbind(pval = cox_model$coefficients[,5])

cox_output <- cbind.data.frame(predictor="cox predictor from the training set",
                               COV=ncol(xtrain),
                               n_features=n_features,
                               pseudo_r2=pseudo_r2,
                               X=rownames(HR),
                               HR)

