library(tidyverse)


### load the data
# methylation data
xtrain <- readRDS("/Cluster_Filespace/Marioni_Group/Yufei/data/methyl_train.rds")
xtest <- readRDS("/Cluster_Filespace/Marioni_Group/Yufei/data/methyl_test.rds")
# phenotype data
w1 <- read.csv("/Cluster_Filespace/Marioni_Group/Yufei/data/W1_Data_Yufei.csv", header=TRUE)
w3 <- read.csv("/Cluster_Filespace/Marioni_Group/Yufei/data/W3_Data_Yufei.csv", header=TRUE)

# % of missing values
## numbers in the training set
length(which(is.na(xtrain))) # 25984
## % in the training set
length(which(is.na(xtrain)))/length(xtrain) # 1.27852e-05
## numbers in the test set
length(which(is.na(xtest))) # 5
## % in the test set
length(which(is.na(xtest)))/length(xtest) #2.820114e-09


# Missing-value imputation (Mean imputation) in x training and testing datasets
for (i in 1:ncol(xtrain)) {
  xtrain[is.na(xtrain[,i]),i] <- mean(xtrain[,i], na.rm=TRUE)
  xtest[is.na(xtest[,i]),i] <- mean(xtest[,i], na.rm=TRUE)
}

# Heart.Disease_Prevalence
w1 <- w1%>%
  mutate(Heart.Disease_Prevalence=replace(Heart.Disease_Incidence,which(is.na(Heart.Disease_Incidence)),1))
w3 <- w3%>%
  mutate(Heart.Disease_Prevalence=replace(Heart.Disease_Incidence,which(is.na(Heart.Disease_Incidence)),1))

# Construct y train
sampletrain_id <- rownames(xtrain)
ytrain_ind <- match(sampletrain_id,w1[,"Sample_Sentrix_ID"])
ytrain_withna <- w1[ytrain_ind,"Heart.Disease_Prevalence"]
agetrain_withna <- w1[ytrain_ind,"age"]
sextrain_withna <- w1[ytrain_ind,"sex"]

# remove ytrain=NA from the training dataset
ytrain <- ytrain_withna[-which(is.na(ytrain_withna))]
xtrain <- xtrain[-which(is.na(ytrain_withna)),]
agetrain <- agetrain_withna[-which(is.na(ytrain_withna))]
sextrain <- sextrain_withna[-which(is.na(ytrain_withna))]

############################################################
# Construct y test
sampletest_id <- rownames(xtest)
ytest_ind <- match(sampletest_id,w3[,"Sample_Sentrix_ID"])
ytest <- w3[ytest_ind,"Heart.Disease_Prevalence"]
agetest <- w3[ytest_ind,"age"]
sextest <- w3[ytest_ind,"sex"]


##########################################################################
# Logistic #
##########################################################################
# create a function to calculate COV for all CpGs and store the CpG ids for the thresholds specified in the n_ftr argument
ftr_sel <- function(x,n_ftr) {
  sd=apply(x,2,sd)
  mean=apply(x,2,mean)
  cov = abs(sd/mean)
  sorted_cov <- sort(cov,index.return=TRUE,decreasing=TRUE)
  CpG_id <- colnames(x)[sorted_cov$ix]
  CpG_id[1:n_ftr]
}
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
    write.csv(active_coef, file = paste0("/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/logistic/coef/", penalty, "coefficient_", outcome, "_", COV[k], ".csv"))
    
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
  write.csv(resultdf,file = paste0("/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/logistic/train_scores/", penalty, "result_", outcome, ".csv"))
  # plot roc (training set)
  ggroc(roctrain,legacy.axes=T)+
    scale_color_hue(name="COV",labels=c("10000","20000","30000","40000","50000","100000","200000","300000","398422"))+
    labs(x="False Postive Rate(1-Specificity)",y="True Positive Rate(Sensitivity)")+
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")
  ggsave(filename=paste0("/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/logistic/train_scores/roctrain",penalty,".png"))
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
  ggsave(PRplot,filename=paste0("/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/logistic/train_scores/PRroctrain",penalty,".png"))
  
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


# set class weights
wts <- rep(1,length(ytrain))
wts[which(ytrain==1)] <- length(ytrain)/(2*length(which(ytrain==1))) #14.0231214
wts[which(ytrain!=1)] <- length(ytrain)/(2*length(which(ytrain!=1))) #0.5184869

# roc, pr curve of best ridge model (test data)
set.seed(123)
fit_ridge <- cv.glmnet(xtrain, ytrain, alpha=0, family="binomial", nfolds=3, weights=wts)
yprob_ridge <- as.numeric(predict(fit_ridge, xtest, type="response", s="lambda.min"))
# roc curve
roctest_ridge <- roc(ytest,yprob_ridge,auc=T)
# precision-recall curve
prtst <- pr.curve(scores.class0=yprob_ridge,weights.class0=ytest,curve=T)
prtest_ridge <- prtst$curve

# roc, pr curve of best lasso model (test data)
ftr <- ftr_sel(xtrain,300000)
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

# roc, pr curve of best enet model (test data)
ftr <- ftr_sel(xtrain,300000)
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

roctest <- list()
roctest[[1]] <- roctest_lasso
roctest[[2]] <- roctest_enet
roctest[[3]] <- roctest_ridge
ggroc(roctest,legacy.axes=T)+
  scale_color_hue(name="Penalty",labels=c("Lasso","enet","Ridge"))+
  labs(x="False Postive Rate(1-Specificity)",y="True Positive Rate(Sensitivity)")+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")
ggsave(filename="/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/logistic/test_scores/roctest.png")
# plot PR curve (training set)
colors=c("Lasso"="blue","enet"="green","Ridge"="red")
PRplot <- ggplot(NULL,aes(x=X1,y=X2))+
  geom_line(data=data.frame(prtest_lasso),aes(color="Lasso"))+
  geom_line(data=data.frame(prtest_enet),aes(color="enet"))+
  geom_line(data=data.frame(prtest_ridge),aes(color="Ridge"))+
  labs(x="Recall",y="Precision")+
  scale_color_manual(name="Penalty",values=colors,breaks=c("Lasso","enet","Ridge"))
ggsave(PRplot,filename="/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/logistic/test_scores/PRroctest.png")


##########################################################################
# Random FOrest #
##########################################################################
library(pROC)
library(PRROC)
library(randomForest)
library(caret)
library(e1071)
# CpG sites selected from Cox models
active_coef1 <- read.csv("/Cluster_Filespace/Marioni_Group/Yufei/output/cox/coef/0.5coef398422_cox.csv")
CpGsites1 <- active_coef1[,1]
# CpG sites selected from logistic classifier
active_coef2 <- read.csv("/Cluster_Filespace/Marioni_Group/Yufei/output/binary/classweights/coef/0.5coefficient_binomial_398422.csv")
CpGsites2 <- active_coef2[-1,][,1]
CpGsites <- union(CpGsites1,CpGsites2)

xtrain_sel <- xtrain[,CpGsites]
xtest_sel <- xtest[,CpGsites]

xtrain_sel <- cbind(xtrain_sel,"age"=agetrain,"sex"=sextrain)
xtest_sel <-cbind(xtest_sel,"age"=agetest,"sex"=sextest)

ytrain_factor = replace(ytrain,which(ytrain==1),"yes")
ytrain_factor = replace(ytrain_factor,which(ytrain_factor==0),"no")
ytest_factor = replace(ytest,which(ytest==1),"yes")
ytest_factor = replace(ytest_factor,which(ytest_factor==0),"no")
seeds <- list()
seeds[[1]] <- rep(432,39)
seeds[[2]] <- rep(432,39)
seeds[[3]] <- rep(432,39)
seeds[[4]] <- 432
control <- trainControl(method='cv',
                        number=3,
                        search='grid',
                        classProbs = TRUE,
                        summaryFunction = multiClassSummary,
                        seeds = seeds)
#create tunegrid with 50 values from 1:50 for mtry to tunning model. Our train function will change number of entry variable at each split according to tunegrid. 
tunegrid <- expand.grid(.mtry = c(3,5,7,9,11,13,15,17,19,21,23,25,27),
                        .splitrule = "gini",
                        .min.node.size = c(1,5,9)) 


set.seed(432)
system.time(
  rf_default <- train(factor(ytrain_factor)~.,
                      data = cbind.data.frame(ytrain_factor,xtrain_sel),
                      method='ranger', 
                      metric='AUC',
                      tuneGrid=tunegrid, 
                      trControl=control)
)
rfresults <- as.data.frame(rf_default$results)
write.csv(rfresults, file = "/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/rf/rfcv.csv")


png(filename = "/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/rf/rfcv.png")
ggplot(rf_default)
dev.off()


yrfprob<-predict(rf_default,xtest_sel,type = "prob")[,2]
rocrf_test <- roc(ytest_factor,yrfprob,plot=T,auc=T)
auc <- auc(rocrf_test)
ggroc(rocrf_test,legacy.axes=T)+
  labs(x="False Postive Rate(1-Specificity)",y="True Positive Rate(Sensitivity)")+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")
ggsave(filename="/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/rf/rf_roctest.png")

# plot PR curve (training set)
pr_rf <- pr.curve(scores.class0 = yrfprob[ytest_factor=="yes"], scores.class1 = yrfprob[ytest_factor=="no"], curve = T)
plot(pr_rf)
prauc <- pr_rf$auc.integral
ggplot(data.frame(pr_rf$curve),aes(x=X1,y=X2,color=X3))+
  geom_line()+
  labs(x="Recall",y="Precision",color="Threshold")
ggsave(filename="/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/rf/rf_prroctest.png")
yrfclass <- predict(rf_default,xtest_sel)
# confusion matrix
res <- list()
res[[1]] <- as.matrix(confusionMatrix(data=yrfclass,reference=factor(ytest_factor),positive="yes"),what = "classes")[1:4,]
res_df <- tibble("method"="random forest", "n features"=ncol(xtrain_sel), "mtry"=rf_default$bestTune$mtry, "AUC"=auc, "PR AUC"=prauc, res)%>%
  unnest_wider(res)
write.csv(res_df, file = "/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/rf/rftestmetrics.csv")

##########################################################################
# SVM(RBF) # 
##########################################################################
# Fit the model 
# Set up Repeated k-fold Cross Validation
train_control <- trainControl(method="cv", number=3, classProbs=T, summaryFunction=multiClassSummary)

set.seed(432)
system.time(
  svmR <- train(xtrain_sel, ytrain_factor, 
                method = "svmRadial", 
                metric="AUC", 
                trControl = train_control, 
                tuneGrid = expand.grid(C = 2^(-5:5),sigma=10^(-3:0)),
                preProcess = c("center","scale"))
)
saveRDS(svmR,file = "/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/svm/svmR.rds")
svmR <- readRDS("/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/svm/svmR.rds")

#View the model
svmRcv <- as.data.frame(svmR$results)
write.csv(svmRcv, file = "/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/svm/svmR_cv.csv")

png(filename = "/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/svm/svmRcv.png")
ggplot(svmR)
dev.off()

svmR <-readRDS("/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/svm/svmR.rds")
# evaluation on the test set
# fitted probability
ysvmrprob <- predict(svmR,xtest_sel,type="prob")$yes
rocsvmR_test <- roc(ytest,ysvmrprob)
ggroc(rocsvmR_test,legacy.axes=T)+
  labs(x="False Postive Rate(1-Specificity)",y="True Positive Rate(Sensitivity)")+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")
ggsave(filename="/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/svm/svmR_roctest.png")
auc_svmR <- auc(rocsvmR_test) 
prtst_svmR <-pr.curve(scores.class0=ysvmrprob[ytest_factor=="yes"],scores.class1=ysvmrprob[ytest_factor=="no"],curve=T)
ggplot(data.frame(prtst_svmR$curve),aes(x=X1,y=X2,color=X3))+
  geom_line()+
  labs(x="Recall",y="Precision",color="Threshold")
ggsave(filename="/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/svm/svmR_prroc.png")
prauc_svmR<-prtst_svmR$auc.integral 

# fitted class
ysvmrclass <- predict(svmR,xtest_sel)
# confusion matrix
res <- list()
res[[1]] <- as.matrix(confusionMatrix(data=ysvmrclass,reference=factor(ytest_factor),positive="yes"),what = "classes")[1:4,]
ressvm_df <- tibble("method"="svm(radial)",
                 "n features"=ncol(xtrain_sel),
                 "sigma"=svmR$bestTune$sigma,
                 "C"=svmR$bestTune$C,
                 "AUC"=auc_svmR, 
                 "PR AUC"=prauc_svmR, res)%>%
  unnest_wider(res)
write.csv(ressvm_df, file = "/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/svm/svmRtestmetrics.csv")

#########################################################################
# svm(linear) #
#########################################################################
# Set up Repeated k-fold Cross Validation
train_control <- trainControl(method="cv", number=3, classProbs=T, summaryFunction=multiClassSummary)

set.seed(432)
system.time(
  svmL <- train(factor(ytrain_factor)~.,
                data = cbind.data.frame(ytrain_factor,xtrain_sel),
                method = "svmLinear", 
                metric="AUC", 
                trControl = train_control, 
                tuneGrid = expand.grid(C = 2^(-5:5)),
                preProcess = c("center","scale"))
)
saveRDS(svmL,file = "/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/svm/svmL.rds")

#View the model
svmLcv <- as.data.frame(svmL$results)
write.csv(svmLcv, file = "/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/svm/svmL_cv.csv")

png(filename = "/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/svm/svmLcv.png")
ggplot(svmL)
dev.off()

svmL <-readRDS("/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/svm/svmL.rds")
# evaluation on the test set
# fitted probability
ysvmLprob <- predict(svmL,xtest_sel,type="prob")$yes
rocsvmL_test <- roc(ytest,ysvmLprob)
ggroc(rocsvmL_test,legacy.axes=T)+
  labs(x="False Postive Rate(1-Specificity)",y="True Positive Rate(Sensitivity)")+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")
ggsave(filename="/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/svm/svmL_roctest.png")
auc_svmL <- auc(rocsvmL_test) 
prtst_svmL <- pr.curve(scores.class0=ysvmLprob[ytest_factor=="yes"],scores.class1=ysvmLprob[ytest_factor=="no"],curve=T)
ggplot(data.frame(prtst_svmL$curve),aes(x=X1,y=X2,color=X3))+
  geom_line()+
  labs(x="Recall",y="Precision",color="Threshold")
ggsave(filename="/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/svm/svmL_prroc.png")
prauc_svmL<-prtst_svmL$auc.integral 

# fitted class
ysvmLclass <- predict(svmL,xtest_sel)
# confusion matrix
res <- list()
res[[1]] <- as.matrix(confusionMatrix(data=ysvmLclass,reference=factor(ytest_factor),positive="yes"),what = "classes")[1:4,]
ressvmL_df <- tibble("method"="svm(Linear)",
                     "n features"=ncol(xtrain_sel),
                     "C"=svmL$bestTune$C,
                     "AUC"=auc_svmL, 
                     "PR AUC"=prauc_svmL, res)%>%
  unnest_wider(res)
write.csv(ressvmL_df, file = "/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/svm/svmLtestmetrics.csv")

##########################################################################
# SVM(Poly) #
##########################################################################
set.seed(432)
system.time(
  svmP <- train(factor(ytrain_factor)~.,
                data = cbind.data.frame(ytrain_factor,xtrain_sel),
                method = "svmPoly", 
                metric="AUC", 
                trControl = train_control, 
                tuneLength = 4,
                preProcess = c("center","scale"))
)
saveRDS(svmP,file = "/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/svm/svmP.rds")

#View the model
svmPcv <- as.data.frame(svmP$results)
write.csv(svmPcv, file = "/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/svm/svmP_cv.csv")

png(filename = "/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/svm/svmPcv.png")
ggplot(svmP)
dev.off()

svmP <- readRDS("/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/svm/svmP.rds")
# evaluation on the test set
# fitted probability
ysvmPprob <- predict(svmP,xtest_sel,type="prob")$yes
rocsvmP_test <- roc(ytest,ysvmPprob)
ggroc(rocsvmP_test,legacy.axes=T)+
  labs(x="False Postive Rate(1-Specificity)",y="True Positive Rate(Sensitivity)")+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")
ggsave(filename="/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/svm/svmP_roctest.png")
auc_svmP <- auc(rocsvmP_test) 
prtst_svmP <- pr.curve(scores.class0=ysvmPprob[ytest_factor=="yes"],scores.class1=ysvmPprob[ytest_factor=="no"],curve=T)
ggplot(data.frame(prtst_svmP$curve),aes(x=X1,y=X2,color=X3))+
  geom_line()+
  labs(x="Recall",y="Precision",color="Threshold")
ggsave(filename="/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/svm/svmP_prroc.png")
prauc_svmP<-prtst_svmP$auc.integral 

# fitted class
ysvmPclass <- predict(svmP,xtest_sel)
# confusion matrix
res <- list()
res[[1]] <- as.matrix(confusionMatrix(data=ysvmPclass,reference=factor(ytest_factor),positive="yes"),what = "classes")[1:4,]
ressvmP_df <- tibble("method"="svm(Poly)",
                     "n features"=ncol(xtrain_sel),
                     "C"=svmP$bestTune$C,
                     "degree"=svmP$bestTune$degree,
                     "scale"=svmP$bestTune$scale,
                     "AUC"=auc_svmP, 
                     "PR AUC"=prauc_svmP, res)%>%
  unnest_wider(res)
write.csv(ressvmP_df, file = "/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/svm/svmPtestmetrics.csv")

##########################################################################
# BART(binary) probit link#
##########################################################################
library(BART)
library(MLmetrics)
## tune hyperparameters
k <- c(1,3,5)
power <- c(0.5,2,4)
ntree <- c(50,100,200)
tune <- list()
i=1
for (n in 1:length(ntree)) {
  for (t in 1:length(k)) {
    for (p in 1:length(power)) {
      auccv <- numeric(3)
      prauccv <- numeric(3)
      set.seed(432)
      idx <- createFolds(factor(ytrain),k=3,list=T)
      
      xt1 <- xtrain_sel[-idx[[1]],]
      yt1 <- ytrain[-idx[[1]]]
      xv1 <- xtrain_sel[idx[[1]],]
      yv1 <- ytrain[idx[[1]]]
      bart_cv1 <- mc.pbart(x.train = xt1,
                           y.train = yt1,
                           k = k[t],
                           power = power[p],
                           ntree = ntree[n],
                           seed = 99,
                           mc.cores = 4)
      ypre1 <- as.numeric(predict(bart_cv1,xv1)$prob.test.mean)
      auccv[1] <- AUC(ypre1,yv1)
      prauccv[1] <- PRAUC(ypre1,yv1)
      
      xt2 <- xtrain_sel[-idx[[2]],]
      yt2 <- ytrain[-idx[[2]]]
      xv2 <- xtrain_sel[idx[[2]],]
      yv2 <- ytrain[idx[[2]]]
      bart_cv2 <- mc.pbart(x.train = xt2,
                           y.train = yt2,
                           k = k[t],
                           power = power[p],
                           ntree = ntree[n],
                           seed = 99,
                           mc.cores = 4)
      ypre2 <- as.numeric(predict(bart_cv2,xv2)$prob.test.mean)
      auccv[2] <- AUC(ypre2,yv2)
      prauccv[2] <- PRAUC(ypre2,yv2)
      
      xt3 <- xtrain_sel[-idx[[3]],]
      yt3 <- ytrain[-idx[[3]]]
      xv3 <- xtrain_sel[idx[[3]],]
      yv3 <- ytrain[idx[[3]]]
      bart_cv3 <- mc.pbart(x.train = xt3,
                           y.train = yt3,
                           k = k[t],
                           power = power[p],
                           ntree = ntree[n],
                           seed = 99,
                           mc.cores = 4)
      ypre3 <- as.numeric(predict(bart_cv3,xv3)$prob.test.mean)
      auccv[3] <- AUC(ypre3,yv3)
      prauccv[3] <- PRAUC(ypre3,yv3)
      
      auc_cv <- mean(auccv)
      prauc_cv <- mean(prauccv)
      tune[[i]] <- c("ntree"=ntree[n],"k"=k[t],"power"=power[p],"AUC"=auc_cv,"PR AUC"=prauc_cv)
      i = i+1
    }
  }
}

tune_table <- tibble(tune)%>%unnest_wider(tune)
write.csv(tune_table,file="/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/bart/bartbinary_cv.csv")
tune_table <- read_csv("/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/bart/bartbinary_cv.csv")
ggcv <- ggplot(tune_table,aes(x=ntree,y=AUC,color=factor(power),shape=factor(k),group=interaction(power,k)))+
  geom_point(size=4)+geom_line()+
  labs(y="AUC(Cross-validation)",color="power",shape="k")
ggsave(ggcv,filename= "/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/bart/bartbinary_cv.png")
filter(.data = tune_table,AUC==max(AUC))
# A tibble: 1 x 5
# ntree   k   power   AUC `PR AUC`
# <dbl> <dbl> <dbl> <dbl>    <dbl>
#  200    3     4   0.816    0.288

bart_binaryprobit <- mc.pbart(x.train = xtrain_sel,
                        y.train = ytrain,
                        k = 3,
                        power = 4,
                        ntree = 200,
                        seed = 99,
                        mc.cores = 4)
yprobitpreprob <- as.numeric(predict(bart_binaryprobit,xtest_sel)$prob.test.mean)
rocbartprobit <- roc(ytest,yprobitpreprob,auc=T)
ggroc(rocbartprobit,legacy.axes = T)+
  labs(x="False Postive Rate(1-Specificity)",y="True Positive Rate(Sensitivity)")+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")
ggsave(filename = "/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/bart/bart_roc.png")
aucbart <- auc(rocbart)
prbartprobit <- pr.curve(scores.class0 = yprobitpreprob,weights.class0 = ytest,curve = T)
ggplot(data.frame(prbart$curve),aes(x=X1,y=X2,color=X3))+
  geom_line()+
  labs(x="Recall",y="Precision",color="Threshold")
ggsave(filename="/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/bart/bart_prroc.png")
praucbart <- prbart$auc.integral

ypreclass <- numeric(length(ypreprob))
ypreclass[ypreprob<0.5] <- 0
ypreclass[ypreprob>=0.5] <- 1
# confusion matrix
res <- list()
res[[1]] <- as.matrix(confusionMatrix(data=factor(ypreclass),reference=factor(ytest),positive="1"),what = "classes")[1:4,]
res_df <- tibble("method"="BART(binary)", "n features"=ncol(xtrain_sel), "k"=3, "power"=4, "ntree"=200, "AUC"=aucbart, "PR AUC"=praucbart, res)%>%
  unnest_wider(res)
write.csv(res_df, file = "/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/bart/bart(binary)testmetrics.csv")

##########################################################################
# BART(binary) logit link#
##########################################################################
library(BART)
library(MLmetrics)
## tune hyperparameters
k <- c(1,3,5)
power <- c(0.5,2,4)
ntree <- c(50,100,200)
tune <- list()
i=1
for (n in 1:length(ntree)) {
  for (t in 1:length(k)) {
    for (p in 1:length(power)) {
      auccv <- numeric(3)
      prauccv <- numeric(3)
      set.seed(432)
      idx <- createFolds(factor(ytrain),k=3,list=T)
      
      xt1 <- xtrain_sel[-idx[[1]],]
      yt1 <- ytrain[-idx[[1]]]
      xv1 <- xtrain_sel[idx[[1]],]
      yv1 <- ytrain[idx[[1]]]
      bart_cv1 <- mc.lbart(x.train = xt1,
                           y.train = yt1,
                           k = k[t],
                           power = power[p],
                           ntree = ntree[n],
                           seed = 99,
                           mc.cores = 4)
      ypre1 <- as.numeric(predict(bart_cv1,xv1)$prob.test.mean)
      auccv[1] <- AUC(ypre1,yv1)
      prauccv[1] <- PRAUC(ypre1,yv1)
      
      xt2 <- xtrain_sel[-idx[[2]],]
      yt2 <- ytrain[-idx[[2]]]
      xv2 <- xtrain_sel[idx[[2]],]
      yv2 <- ytrain[idx[[2]]]
      bart_cv2 <- mc.lbart(x.train = xt2,
                           y.train = yt2,
                           k = k[t],
                           power = power[p],
                           ntree = ntree[n],
                           seed = 99,
                           mc.cores = 4)
      ypre2 <- as.numeric(predict(bart_cv2,xv2)$prob.test.mean)
      auccv[2] <- AUC(ypre2,yv2)
      prauccv[2] <- PRAUC(ypre2,yv2)
      
      xt3 <- xtrain_sel[-idx[[3]],]
      yt3 <- ytrain[-idx[[3]]]
      xv3 <- xtrain_sel[idx[[3]],]
      yv3 <- ytrain[idx[[3]]]
      bart_cv3 <- mc.lbart(x.train = xt3,
                           y.train = yt3,
                           k = k[t],
                           power = power[p],
                           ntree = ntree[n],
                           seed = 99,
                           mc.cores = 4)
      ypre3 <- as.numeric(predict(bart_cv3,xv3)$prob.test.mean)
      auccv[3] <- AUC(ypre3,yv3)
      prauccv[3] <- PRAUC(ypre3,yv3)
      
      auc_cv <- mean(auccv)
      prauc_cv <- mean(prauccv)
      tune[[i]] <- c("ntree"=ntree[n],"k"=k[t],"power"=power[p],"AUC"=auc_cv,"PR AUC"=prauc_cv)
      i = i+1
    }
  }
}

tune_table <- tibble(tune)%>%unnest_wider(tune)
write.csv(tune_table,file="/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/bart/bartbinarylogit_cv.csv")
tune_table <- read_csv("/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/bart/bartbinarylogit_cv.csv")
ggcv <- ggplot(tune_table,aes(x=ntree,y=AUC,color=factor(power),shape=factor(k),group=interaction(power,k)))+
  geom_point(size=4)+geom_line()+
  labs(y="AUC(Cross-validation)",color="power",shape="k")
ggsave(ggcv,filename= "/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/bart/bartbinarylogit_cv.png")
filter(.data = tune_table,AUC==max(AUC))
# A tibble: 1 x 5
# ntree  k    power   AUC `PR AUC`
# <dbl> <dbl> <dbl> <dbl>    <dbl>
# 200     3     2   0.817    0.284

bart_binary <- mc.lbart(x.train = xtrain_sel,
                        y.train = ytrain,
                        k = 3,
                        power = 2,
                        ntree = 200,
                        seed = 99,
                        mc.cores = 4)
ypreprob <- as.numeric(predict(bart_binary,xtest_sel)$prob.test.mean)
rocbart <- roc(ytest,ypreprob,auc=T)
ggroc(rocbart,legacy.axes = T)+
  labs(x="False Postive Rate(1-Specificity)",y="True Positive Rate(Sensitivity)")+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")
ggsave(filename = "/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/bart/bartlogit_roc.png")
aucbart <- auc(rocbart)
prbart <- pr.curve(scores.class0 = ypreprob,weights.class0 = ytest,curve = T)
ggplot(data.frame(prbart$curve),aes(x=X1,y=X2,color=X3))+
  geom_line()+
  labs(x="Recall",y="Precision",color="Threshold")
ggsave(filename="/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/bart/bartlogit_prroc.png")
praucbart <- prbart$auc.integral

ypreclass <- numeric(length(ypreprob))
ypreclass[ypreprob<0.5] <- 0
ypreclass[ypreprob>=0.5] <- 1
# confusion matrix
res <- list()
res[[1]] <- as.matrix(confusionMatrix(data=factor(ypreclass),reference=factor(ytest),positive="1"),what = "classes")[1:4,]
res_df <- tibble("method"="logit BART(binary)", "n features"=ncol(xtrain_sel), "k"=3, "power"=2, "ntree"=200, "AUC"=aucbart, "PR AUC"=praucbart, res)%>%
  unnest_wider(res)
write.csv(res_df, file = "/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/bart/bart(binarylogit)testmetrics.csv")

roclist <- list()
roclist[[1]] <- rocrf_test
roclist[[2]] <- rocsvmL_test
roclist[[3]] <- rocsvmR_test
roclist[[4]] <- rocsvmP_test
roclist[[5]] <- rocbartprobit
roclist[[6]] <- rocbart #logit bart

ggroc(roclist,legacy.axes = T)+
  scale_color_hue(name="Models",labels=c("RF","SVM(Linear)","SVM(RBF)","SVM(Poly)","probit BART","logit BART"))+
  labs(x="False Postive Rate(1-Specificity)",y="True Positive Rate(Sensitivity)")+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")
ggsave(filename = "/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/roctestML.png")

colors=c("RF"="black","SVM(Linear)"="brown","SVM(RBF)"="purple","SVM(Poly)"="blue","probit BART"="green","logit BART"="yellow")
PRplot <- ggplot(NULL,aes(x=X1,y=X2))+
  geom_line(data=data.frame(pr_rf$curve),aes(color="RF"))+
  geom_line(data=data.frame(prtst_svmL$curve),aes(color="SVM(Linear)"))+
  geom_line(data=data.frame(prtst_svmR$curve),aes(color="SVM(RBF)"))+
  geom_line(data=data.frame(prtst_svmP$curve),aes(color="SVM(Poly)"))+
  geom_line(data=data.frame(prbartprobit$curve),aes(color="probit BART"))+
  geom_line(data=data.frame(prbart$curve),aes(color="logit BART"))+
  labs(x="Recall",y="Precision")+
  scale_color_manual(name="COV",values=colors,breaks=c("RF","SVM(Linear)","SVM(RBF)","SVM(Poly)","probit BART","logit BART"))
ggsave(PRplot,filename="/Cluster_Filespace/Marioni_Group/Yufei/output/prevalence/prroctestML.png")
