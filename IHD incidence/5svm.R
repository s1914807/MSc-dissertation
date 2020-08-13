library(caret)
library(e1071)
library(pROC)
library(PRROC)
library(tidyverse)
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
# Set up Repeated k-fold Cross Validation
train_control <- trainControl(method="cv", number=3, classProbs=T, summaryFunction=multiClassSummary)

#########################################################################
# svm(RBF) #
#########################################################################
# Fit the model 
set.seed(432)
system.time(
  svmR <- train(xtrain_sel, ytrain_factor, 
                method = "svmRadial", 
                metric="AUC", 
                trControl = train_control, 
                tuneGrid = expand.grid(C = 2^(-5:5),sigma=10^(-3:0)),
                preProcess = c("center","scale"))
)
saveRDS(svmR,file = "/Cluster_Filespace/Marioni_Group/Yufei/output/svm/svmR.rds")

#View the model
svmRcv <- as.data.frame(svmR$results)
write.csv(svmRcv, file = "/Cluster_Filespace/Marioni_Group/Yufei/output/svm/svmR_cv.csv")

png(filename = "/Cluster_Filespace/Marioni_Group/Yufei/output/svm/svmRcv.png")
ggplot(svmR)
dev.off()

svmR <-readRDS("/Cluster_Filespace/Marioni_Group/Yufei/output/svm/svmR.rds")
# evaluation on the test set
# fitted probability
ysvmrprob <- predict(svmR,xtest_sel,type="prob")$yes
rocsvmR_test <- roc(ytest,ysvmrprob)
ggroc(rocsvmR_test,legacy.axes=T)+
  labs(x="False Postive Rate(1-Specificity)",y="True Positive Rate(Sensitivity)")+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")
ggsave(filename="/Cluster_Filespace/Marioni_Group/Yufei/output/svm/svmR_roctest.png")

auc_svmR <- auc(rocsvmR_test) 
prtst_svmR <- pr.curve(scores.class0=ysvmrprob[ytest_factor=="yes"],scores.class1=ysvmrprob[ytest_factor=="no"],curve=T)
ggplot(data.frame(prtst_svmR$curve),aes(x=X1,y=X2,color=X3))+
  geom_line()+
  labs(x="Recall",y="Precision",color="Threshold")
ggsave(filename="/Cluster_Filespace/Marioni_Group/Yufei/output/svm/svmR_prroc.png")

prauc_svmR<-prtst_svmR$auc.integral 

# fitted class
ysvmrclass <- predict(svmR,xtest_sel)
# confusion matrix
res <- list()
res[[1]] <- as.matrix(confusionMatrix(data=ysvmrclass,reference=factor(ytest_factor),positive="yes"),what = "classes")[1:4,]
res_df <- tibble("method"="svm(radial)",
                 "n features"=ncol(xtrain_sel),
                 "sigma"=svmR$bestTune$sigma,
                 "C"=svmR$bestTune$C,
                 "AUC"=auc_svmR, 
                 "PR AUC"=prauc_svmR, res)%>%
  unnest_wider(res)
write.csv(res_df, file = "/Cluster_Filespace/Marioni_Group/Yufei/output/svm/svmRtestmetrics.csv")


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
saveRDS(svmL,file = "/Cluster_Filespace/Marioni_Group/Yufei/output/svm/svmL.rds")

#View the model
svmLcv <- as.data.frame(svmL$results)
write.csv(svmLcv, file = "/Cluster_Filespace/Marioni_Group/Yufei/output/svm/svmL_cv.csv")

png(filename = "/Cluster_Filespace/Marioni_Group/Yufei/output/svm/svmLcv.png")
ggplot(svmL)
dev.off()

svmL <- readRDS("/Cluster_Filespace/Marioni_Group/Yufei/output/svm/svmL.rds")
# evaluation on the test set
# fitted probability
ysvmLprob <- predict(svmL,xtest_sel,type="prob")$yes
rocsvmL_test <- roc(ytest,ysvmLprob)
ggroc(rocsvmL_test,legacy.axes=T)+
  labs(x="False Postive Rate(1-Specificity)",y="True Positive Rate(Sensitivity)")+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")
ggsave(filename="/Cluster_Filespace/Marioni_Group/Yufei/output/svm/svmL_roctest.png")
auc_svmL <- auc(rocsvmL_test) 
prtst_svmL<-pr.curve(scores.class0=ysvmLprob[ytest_factor=="yes"],scores.class1=ysvmLprob[ytest_factor=="no"],curve=T)
ggplot(data.frame(prtst_svmL$curve),aes(x=X1,y=X2,color=X3))+
  geom_line()+
  labs(x="Recall",y="Precision",color="Threshold")
ggsave(filename="/Cluster_Filespace/Marioni_Group/Yufei/output/svm/svmL_prroc.png")
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
write.csv(ressvmL_df, file = "/Cluster_Filespace/Marioni_Group/Yufei/output/svm/svmLtestmetrics.csv")

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
saveRDS(svmP,file = "/Cluster_Filespace/Marioni_Group/Yufei/output/svm/svmP.rds")

#View the model
svmPcv <- as.data.frame(svmP$results)
write.csv(svmPcv, file = "/Cluster_Filespace/Marioni_Group/Yufei/output/svm/svmP_cv.csv")

png(filename = "/Cluster_Filespace/Marioni_Group/Yufei/output/svm/svmPcv.png")
ggplot(svmP)
dev.off()

svmP <-readRDS("/Cluster_Filespace/Marioni_Group/Yufei/output/svm/svmP.rds")
# evaluation on the test set
# fitted probability
ysvmPprob <- predict(svmP,xtest_sel,type="prob")$yes
rocsvmP_test <- roc(ytest,ysvmPprob)
ggroc(rocsvmP_test,legacy.axes=T)+
  labs(x="False Postive Rate(1-Specificity)",y="True Positive Rate(Sensitivity)")+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")
ggsave(filename="/Cluster_Filespace/Marioni_Group/Yufei/output/svm/svmP_roctest.png")
auc_svmP <- auc(rocsvmP_test) 
prtst_svmP <- pr.curve(scores.class0=ysvmPprob[ytest_factor=="yes"],scores.class1=ysvmPprob[ytest_factor=="no"],curve=T)
ggplot(data.frame(prtst_svmP$curve),aes(x=X1,y=X2,color=X3))+
  geom_line()+
  labs(x="Recall",y="Precision",color="Threshold")
ggsave(filename="/Cluster_Filespace/Marioni_Group/Yufei/output/svm/svmP_prroc.png")
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
write.csv(ressvmP_df, file = "/Cluster_Filespace/Marioni_Group/Yufei/output/svm/svmPtestmetrics.csv")
