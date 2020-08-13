library(tidyverse)
library(ggplot2)
library(pROC)
library(PRROC)
library(caret)
library(BART)
library(MLmetrics)
# CpG sites selected from Cox models
active_coef1 <- read.csv("/Cluster_Filespace/Marioni_Group/Yufei/output/cox/coef/0.5coef398422_cox.csv")
CpGsites1 <- active_coef1[,1]
# CpG sites selected from logistic classifier
active_coef2 <- read.csv("/Cluster_Filespace/Marioni_Group/Yufei/output/binary/classweights/coef/0.5coefficient_binomial_398422.csv")
CpGsites2 <- active_coef2[-1,][,1]
CpGsites <- union(CpGsites1,CpGsites2)
xtrain_sel <- xtrain[,CpGsites]%>%cbind("age"=agetrain,"sex"=sextrain)
xtest_sel <- xtest[,CpGsites]%>%cbind("age"=agetest,"sex"=sextest)


# probit BART
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
write.csv(tune_table,file="/Cluster_Filespace/Marioni_Group/Yufei/output/bart/bartbinary_cv.csv")
tune_table <- read_csv("/Cluster_Filespace/Marioni_Group/Yufei/output/bart/bartbinary_cv.csv")
ggcv <- ggplot(tune_table,aes(x=ntree,y=AUC,color=factor(power),shape=factor(k),group=interaction(power,k)))+
  geom_point(size=4)+geom_line()+
  labs(y="AUC(Cross-validation)",color="power",shape="k")
ggsave(ggcv,filename= "/Cluster_Filespace/Marioni_Group/Yufei/output/bart/bartbinary_cv.png")
filter(.data = tune_table,AUC==max(AUC))
# A tibble: 1 x 5
# ntree     k power   AUC `PR AUC`
# <dbl> <dbl> <dbl> <dbl>    <dbl>
# 200     3   0.5   0.828    0.166

bart_binary <- mc.pbart(x.train = xtrain_sel,
                        y.train = ytrain,
                        k = 3,
                        power = 0.5,
                        ntree = 200,
                        seed = 99,
                        mc.cores = 4)
ypreprob <- as.numeric(predict(bart_binary,xtest_sel)$prob.test.mean)
rocbart <- roc(ytest,ypreprob,auc=T)
ggroc(rocbart,legacy.axes = T)+
  labs(x="False Postive Rate(1-Specificity)",y="True Positive Rate(Sensitivity)")+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")
ggsave(filename = "/Cluster_Filespace/Marioni_Group/Yufei/output/bart/bart_roc.png")
aucbart <- auc(rocbart)
prbart <- pr.curve(scores.class0 = ypreprob,weights.class0 = ytest,curve = T)
ggplot(data.frame(prbart$curve),aes(x=X1,y=X2,color=X3))+
  geom_line()+
  labs(x="Recall",y="Precision",color="Threshold")
ggsave(filename="/Cluster_Filespace/Marioni_Group/Yufei/output/bart/bart_prroc.png")
praucbart <- prbart$auc.integral

ypreclass <- numeric(length(ypreprob))
ypreclass[ypreprob<0.5] <- 0
ypreclass[ypreprob>=0.5] <- 1
# confusion matrix
res <- list()
res[[1]] <- as.matrix(confusionMatrix(data=factor(ypreclass),reference=factor(ytest),positive="1"),what = "classes")[1:4,]
res_df <- tibble("method"="BART(binary)", "n features"=ncol(xtrain_sel), "k"=3, "power"=0.5, "ntree"=200, "AUC"=aucbart, "PR AUC"=praucbart, res)%>%
  unnest_wider(res)
write.csv(res_df, file = "/Cluster_Filespace/Marioni_Group/Yufei/output/bart/bart(binary)testmetrics.csv")


##logit BART
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
write.csv(tune_table,file="/Cluster_Filespace/Marioni_Group/Yufei/output/bart/logit/bartbinary_cv.csv")
tune_table <- read_csv("/Cluster_Filespace/Marioni_Group/Yufei/output/bart/logit/bartbinary_cv.csv")
ggcv <- ggplot(tune_table,aes(x=ntree,y=AUC,color=factor(power),shape=factor(k),group=interaction(power,k)))+
  geom_point(size=4)+geom_line()+
  labs(y="AUC(Cross-validation)",color="power",shape="k")
ggsave(ggcv,filename= "/Cluster_Filespace/Marioni_Group/Yufei/output/bart/logit/bartbinary_cv.png")
filter(.data = tune_table,AUC==max(AUC))
# A tibble: 1 x 5
# ntree     k power   AUC `PR AUC`
# <dbl> <dbl> <dbl> <dbl>    <dbl>
# 200     3     2   0.830    0.171

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
ggsave(filename = "/Cluster_Filespace/Marioni_Group/Yufei/output/bart/logit/bart_roc.png")
aucbart <- auc(rocbart)
prbart <- pr.curve(scores.class0 = ypreprob,weights.class0 = ytest,curve = T)
ggplot(data.frame(prbart$curve),aes(x=X1,y=X2,color=X3))+
  geom_line()+
  labs(x="Recall",y="Precision",color="Threshold")
ggsave(filename="/Cluster_Filespace/Marioni_Group/Yufei/output/bart/logit/bart_prroc.png")
praucbart <- prbart$auc.integral

ypreclass <- numeric(length(ypreprob))
ypreclass[ypreprob<0.5] <- 0
ypreclass[ypreprob>=0.5] <- 1
# confusion matrix
res <- list()
res[[1]] <- as.matrix(confusionMatrix(data=factor(ypreclass),reference=factor(ytest),positive="1"),what = "classes")[1:4,]
res_df <- tibble("method"="BART(binary)", "n features"=ncol(xtrain_sel), "k"=3, "power"=0.5, "ntree"=200, "AUC"=aucbart, "PR AUC"=praucbart, res)%>%
  unnest_wider(res)
write.csv(res_df, file = "/Cluster_Filespace/Marioni_Group/Yufei/output/bart/logit/bart(binary)testmetrics.csv")


