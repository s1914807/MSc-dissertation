library(tidyverse)
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
write.csv(rfresults, file = "/Cluster_Filespace/Marioni_Group/Yufei/output/rf/rfcv.csv")


png(filename = "/Cluster_Filespace/Marioni_Group/Yufei/output/rf/rfcv.png")
ggplot(rf_default)
dev.off()


yrfprob<-predict(rf_default,xtest_sel,type = "prob")[,2]
rocrf_test <- roc(ytest_factor,yrfprob,plot=T,auc=T)
auc <- auc(rocrf_test)
ggroc(rocrf_test,legacy.axes=T)+
  labs(x="False Postive Rate(1-Specificity)",y="True Positive Rate(Sensitivity)")+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")
ggsave(filename="/Cluster_Filespace/Marioni_Group/Yufei/output/rf/rf_roctest.png")

# plot PR curve 
pr_rf <- pr.curve(scores.class0 = yrfprob[ytest_factor=="yes"], scores.class1 = yrfprob[ytest_factor=="no"], curve = T)
plot(pr_rf)
prauc <- pr_rf$auc.integral
ggplot(data.frame(pr_rf$curve),aes(x=X1,y=X2,color=X3))+
  geom_line()+
  labs(x="Recall",y="Precision",color="Threshold")
ggsave(filename="/Cluster_Filespace/Marioni_Group/Yufei/output/rf/rf_prroctest.png")
yrfclass <- predict(rf_default,xtest_sel)
# confusion matrix
res <- list()
res[[1]] <- as.matrix(confusionMatrix(data=yrfclass,reference=factor(ytest_factor),positive="yes"),what = "classes")[1:4,]
res_df <- tibble("method"="random forest", "n features"=ncol(xtrain_sel), "mtry"=rf_default$bestTune$mtry, "AUC"=auc, "PR AUC"=prauc, res)%>%
  unnest_wider(res)
write.csv(res_df, file = "/Cluster_Filespace/Marioni_Group/Yufei/output/rf/rftestmetrics.csv")
