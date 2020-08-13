# data preprocessing

### load the data
# methylation data
xtrain <- readRDS("/Cluster_Filespace/Marioni_Group/Yufei/data/methyl_train.rds")
xtest <- readRDS("/Cluster_Filespace/Marioni_Group/Yufei/data/methyl_test.rds")
# phenotype data
w1 <- read.csv("/Cluster_Filespace/Marioni_Group/Yufei/data/W1_Data_Yufei.csv", header=TRUE)
w3 <- read.csv("/Cluster_Filespace/Marioni_Group/Yufei/data/W3_Data_Yufei.csv", header=TRUE)


# Missing-value imputation (Mean imputation) in x training and testing datasets
for (i in 1:ncol(xtrain)) {
  xtrain[is.na(xtrain[,i]),i] <- mean(xtrain[,i], na.rm=TRUE)
  xtest[is.na(xtest[,i]),i] <- mean(xtest[,i], na.rm=TRUE)
}

# for people whose "Heart.Disease_Age_Event"<"age", lablel Heart.Disease_Incidence=NA
w1 <- w1%>%
  mutate(Heart.Disease_Incidence=replace(Heart.Disease_Incidence,which(Heart.Disease_Age_Event<age&Heart.Disease_Incidence==1),NA))
w3 <- w3%>%
  mutate(Heart.Disease_Incidence=replace(Heart.Disease_Incidence,which(Heart.Disease_Age_Event<age&Heart.Disease_Incidence==1),NA))

# Construct y train
sampletrain_id <- rownames(xtrain)
ytrain_ind <- match(sampletrain_id,w1[,"Sample_Sentrix_ID"])
ytrain_withna <- w1[ytrain_ind,"Heart.Disease_Incidence"]
agetrain_withna <- w1[ytrain_ind,"age"]
sextrain_withna <- w1[ytrain_ind,"sex"]
bmitrain_withna <- w1[ytrain_ind,"bmi"]

# construct age_event_train
for (i in 1:nrow(w1)){
  if(is.na(w1[i,"Heart.Disease_Age_Event"])){
    w1[i,"Heart.Disease_Age_Event"] <- w1[i,"aged"]
  }
}
age_event_train_withna <- w1[ytrain_ind,"Heart.Disease_Age_Event"]-w1[ytrain_ind,"age"]

# remove individuals with Heart.Disease_Incidence=NA from the training dataset
ytrain <- ytrain_withna[-which(is.na(ytrain_withna))]
xtrain <- xtrain[-which(is.na(ytrain_withna)),]
agetrain <- agetrain_withna[-which(is.na(ytrain_withna))]
sextrain <- sextrain_withna[-which(is.na(ytrain_withna))]
age_event_train_withzero <- age_event_train_withna[-which(is.na(ytrain_withna))]
bmitrain <- bmitrain_withna[-which(is.na(ytrain_withna))]
############################################################


# Construct y test
sampletest_id <- rownames(xtest)
ytest_ind <- match(sampletest_id,w3[,"Sample_Sentrix_ID"])
ytest_withna <- w3[ytest_ind,"Heart.Disease_Incidence"]
agetest_withna <- w3[ytest_ind,"age"]
sextest_withna <- w3[ytest_ind,"sex"]
bmitest_withna <- w3[ytest_ind,"bmi"]

for (i in 1:nrow(w3)){
  if(is.na(w3[i,"Heart.Disease_Age_Event"])){
    w3[i,"Heart.Disease_Age_Event"] <- w3[i,"aged"]
  }
}
age_event_test_withna <- w3[ytest_ind,"Heart.Disease_Age_Event"]-w3[ytest_ind,"age"]

# remove individuals with Heart.Disease_Incidence=NA from the testing dataset
ytest <- ytest_withna[-which(is.na(ytest_withna))]
xtest <- xtest[-which(is.na(ytest_withna)),]
agetest <- agetest_withna[-which(is.na(ytest_withna))]
sextest <- sextest_withna[-which(is.na(ytest_withna))]
age_event_test_withneg <- age_event_test_withna[-which(is.na(ytest_withna))]
bmitest <- bmitest_withna[-which(is.na(ytest_withna))]
####################################################

#remove age_event_train_withzero<=0
age_event_train <- age_event_train_withzero[-which(age_event_train_withzero<=0)]
ytrain <- ytrain[-which(age_event_train_withzero<=0)]
xtrain <- xtrain[-which(age_event_train_withzero<=0),]
agetrain <- agetrain[-which(age_event_train_withzero<=0)]
sextrain <- sextrain[-which(age_event_train_withzero<=0)]
bmitrain <- bmitrain[-which(age_event_train_withzero<=0)]

# remove negative time in test data
age_event_test <- age_event_test_withneg[which(age_event_test_withneg>0)]
ytest <- ytest[which(age_event_test_withneg>0)]
xtest <- xtest[which(age_event_test_withneg>0),]
agetest <- agetest[which(age_event_test_withneg>0)]
sextest <- sextest[which(age_event_test_withneg>0)]
bmitest <- bmitest[which(age_event_test_withneg>0)]
###########################################################################
