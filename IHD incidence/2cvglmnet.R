library(pROC)
library(PRROC)
library(ggplot2)
library(MLmetrics)
library(caret)
library(glmnet)


pip_logistic <- function(xtrain,ytrain,xtest,ytest,penalty,outcome,nfolds_cv,COV) {
  auc <- numeric(length(COV))
  prauc <- numeric(length(COV))
  result <- list()
  roctrain <- list()
  prtrain <- list()
  
  
  auc_test <- numeric(length(COV))
  prauc_test <- numeric(length(COV))
  res <- list()
  roctest<- list()
  prtest <- list()
  
  # set class weights
  wts <- rep(1,length(ytrain))
  wts[which(ytrain==1)] <- length(ytrain)/(2*length(which(ytrain==1))) #14.0231214
  wts[which(ytrain!=1)] <- length(ytrain)/(2*length(which(ytrain!=1))) #0.5184869
  
  
  for (k in 1:length(COV)) {
    ftr <- ftr_sel(xtrain,COV[k])
    # Construct X training subset
    Xtrain_sub <- xtrain[,ftr]
    # fit the model
    set.seed(123)
    fit <- cv.glmnet(Xtrain_sub, ytrain, alpha=penalty, family=outcome, nfolds=nfolds_cv, weights=wts)
    # extract coefficients
    coef <- coef(fit, s="lambda.min")
    active_coef <- as.data.frame(coef[which(coef!=0),])
    write.csv(active_coef, file = paste0("/Cluster_Filespace/Marioni_Group/Yufei/output/binary/classweights/coef/", penalty, "coefficient_", outcome, "_", COV[k], ".csv"))
    
    # results of the training set
    # predicted class for the training set
    yfittedclass <- as.numeric(predict(fit, Xtrain_sub, type="class", s="lambda.min"))
    # sensitivity, specificity, ppv and npv from confusion matrix
    result[[k]] <- as.matrix(confusionMatrix(as.factor(yfittedclass),as.factor(ytrain),positive="1"),what = "classes")[1:4,]
    # fitted probabilities for the training set
    yfittedprob <- as.numeric(predict(fit, Xtrain_sub, type="response", s="lambda.min"))
    # roc curve
    roctrain[[k]] <- roc(ytrain,yfittedprob,auc=T)
    # auc
    auc[k] <- auc(roctrain[[k]])[1]
    # precision-recall curve
    pr <- pr.curve(scores.class0=yfittedprob,weights.class0=ytrain,curve=T)
    prtrain[[k]] <- pr$curve
    # PR auc
    prauc[k] <- pr$auc.integral
    
    # results of the test set
    # Construct X testing subset
    Xtest_sub <- xtest[,ftr]
    # predict the class for the test set
    ypre <- as.numeric(predict(fit, Xtest_sub, type="class", s="lambda.min"))
    # sensitivity, specificity, ppv and npv from confusion matrix
    res[[k]] <- as.matrix(confusionMatrix(as.factor(ypre),as.factor(ytest),positive="1"),what = "classes")[1:4,]
    # fitted probabilities for the test set
    yprob <- as.numeric(predict(fit, Xtest_sub, type="response", s="lambda.min"))
    # roc curve
    roctest[[k]] <- roc(ytest,yprob,auc=T)
    # auc
    auc_test[k] <- auc(roctest[[k]])[1]
    # precision-recall curve
    prtst <- pr.curve(scores.class0=yprob,weights.class0=ytest,curve=T)
    prtest[[k]] <- prtst$curve
    # PR auc
    prauc_test[k] <- prtst$auc.integral
  }
  
  # results and best model
  m <- which(auc==max(auc))[1]
  cat("The best model: COV=",COV[m])
  resultdf <- tibble("COV"=COV,"AUC"=auc,"PR AUC"=prauc,result)%>%
    unnest_wider(result)
  write.csv(resultdf,file = paste0("/Cluster_Filespace/Marioni_Group/Yufei/output/binary/classweights/train_scores/", penalty, "result_", outcome, ".csv"))
  # plot roc (training set)
  ggroc(roctrain,legacy.axes=T)+
    scale_color_hue(name="COV",labels=c("10000","20000","30000","40000","50000","100000","200000","300000","398422"))+
    labs(x="False Postive Rate(1-Specificity)",y="True Positive Rate(Sensitivity)")+
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")
  ggsave(filename=paste0("/Cluster_Filespace/Marioni_Group/Yufei/output/binary/classweights/train_scores/roctrain",penalty,".png"))
  # plot PR curve (training set)
  colors=c("10000"="black","20000"="brown","30000"="purple","40000"="blue","50000"="green","100000"="yellow","200000"="orange","300000"="red","398422"="pink")
  PRplot <- ggplot(NULL,aes(x=X1,y=X2))+
    geom_line(data=data.frame(prtrain[[1]]),aes(color="10000"))+
    geom_line(data=data.frame(prtrain[[2]]),aes(color="20000"))+
    geom_line(data=data.frame(prtrain[[3]]),aes(color="30000"))+
    geom_line(data=data.frame(prtrain[[4]]),aes(color="40000"))+
    geom_line(data=data.frame(prtrain[[5]]),aes(color="50000"))+
    geom_line(data=data.frame(prtrain[[6]]),aes(color="100000"))+
    geom_line(data=data.frame(prtrain[[7]]),aes(color="200000"))+
    geom_line(data=data.frame(prtrain[[8]]),aes(color="300000"))+
    geom_line(data=data.frame(prtrain[[9]]),aes(color="398422"))+
    labs(x="Recall",y="Precision")+
    scale_color_manual(name="COV",values=colors,breaks=c("10000","20000","30000","40000","50000","100000","200000","300000","398422"))
  ggsave(PRplot,filename=paste0("/Cluster_Filespace/Marioni_Group/Yufei/output/binary/classweights/train_scores/PRroctrain",penalty,".png"))
    
  
  
  # evaluation results with test data
  eval <- tibble("COV"=COV,"AUC"=auc_test,"PR AUC"=prauc_test,res)%>%
    unnest_wider(res)
  write.csv(eval,file = paste0("/Cluster_Filespace/Marioni_Group/Yufei/output/binary/classweights/test_scores/", penalty, "eval_",outcome,"_", COV[m], ".csv"))
  # roctestbest <- roctest[[m]]
  # return(roctestbest)
  # ggroc(roctestbest)+
  #   scale_color_hue(name="COV",labels=COV[m])
  # ggsave(filename=paste0("/Cluster_Filespace/Marioni_Group/Yufei/output/binary/classweights/test_scores/roctest",penalty,".png"))
  # 
}



COV <- c(10000,20000,30000,40000,50000,100000,200000,300000,ncol(xtrain))

system.time(
  pip_logistic(xtrain,ytrain,xtest,ytest,penalty=1,outcome="binomial",nfolds_cv=3,COV)
)
gc()
system.time(
  pip_logistic(xtrain,ytrain,xtest,ytest,penalty=0.5,outcome="binomial",nfolds_cv=3,COV)
)
gc()
system.time(
  pip_logistic(xtrain,ytrain,xtest,ytest,penalty=0,outcome="binomial",nfolds_cv=3,COV)
)
gc()
##################################################################

##################################################################
# set class weights
wts <- rep(1,length(ytrain))
wts[which(ytrain==1)] <- length(ytrain)/(2*length(which(ytrain==1))) #14.0231214
wts[which(ytrain!=1)] <- length(ytrain)/(2*length(which(ytrain!=1))) #0.5184869

# roc, pr curve of best enet model (test data)
ftr <- ftr_sel(xtrain,ncol(xtrain))
# Construct X training subset
Xtrain_sub <- xtrain[,ftr]
set.seed(123)
fit_enet <- cv.glmnet(Xtrain_sub, ytrain, alpha=0.5, family="binomial", nfolds=3, weights=wts)
yprob_enet <- as.numeric(predict(fit_enet, xtest[,ftr], type="response", s="lambda.min"))
# roc curve
roctest_enet <- roc(ytest,yprob_enet,auc=T)
# precision-recall curve
prtst <- pr.curve(scores.class0=yprob_enet,weights.class0=ytest,curve=T)
prtest_enet <- prtst$curve

# roc, pr curve of best ridge model (test data)
ftr <- ftr_sel(xtrain,300000)
# Construct X training subset
Xtrain_sub <- xtrain[,ftr]
set.seed(123)
fit_ridge <- cv.glmnet(Xtrain_sub, ytrain, alpha=0, family="binomial", nfolds=3, weights=wts)
yprob_ridge <- as.numeric(predict(fit_ridge, xtest[,ftr], type="response", s="lambda.min"))
# roc curve
roctest_ridge <- roc(ytest,yprob_ridge,auc=T)
# precision-recall curve
prtst <- pr.curve(scores.class0=yprob_ridge,weights.class0=ytest,curve=T)
prtest_ridge <- prtst$curve

# roc, pr curve of best lasso model (test data)
ftr <- ftr_sel(xtrain,100000)
# Construct X training subset
Xtrain_sub <- xtrain[,ftr]
set.seed(123)
fit_lasso <- cv.glmnet(Xtrain_sub, ytrain, alpha=1, family="binomial", nfolds=3, weights=wts)
yprob_lasso <- as.numeric(predict(fit_lasso, xtest[,ftr], type="response", s="lambda.min"))
# roc curve
roctest_lasso <- roc(ytest,yprob_lasso,auc=T)
# precision-recall curve
prtst <- pr.curve(scores.class0=yprob_lasso,weights.class0=ytest,curve=T)
prtest_lasso <- prtst$curve



roctest <- list()
roctest[[1]] <- roctest_lasso
roctest[[2]] <- roctest_enet
roctest[[3]] <- roctest_ridge
ggroc(roctest,legacy.axes=T)+
  scale_color_hue(name="Penalty",labels=c("Lasso","enet","Ridge"))+
  labs(x="False Postive Rate(1-Specificity)",y="True Positive Rate(Sensitivity)")+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")
ggsave(filename="/Cluster_Filespace/Marioni_Group/Yufei/output/binary/classweights/test_scores/roctest.png")
# plot PR curve (training set)
colors=c("Lasso"="blue","enet"="green","Ridge"="red")
PRplot <- ggplot(NULL,aes(x=X1,y=X2))+
  geom_line(data=data.frame(prtest_lasso),aes(color="Lasso"))+
  geom_line(data=data.frame(prtest_enet),aes(color="enet"))+
  geom_line(data=data.frame(prtest_ridge),aes(color="Ridge"))+
  labs(x="Recall",y="Precision")+
  scale_color_manual(name="Penalty",values=colors,breaks=c("Lasso","enet","Ridge"))
ggsave(PRplot,filename="/Cluster_Filespace/Marioni_Group/Yufei/output/binary/classweights/test_scores/PRroctest.png")
