library(dplyr)

anno <- readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds")
probes <- rownames(anno)[anno$Methyl450_Loci!=""]

# convert m-value to beta-value 
m2beta <- function(m) {
  beta <- 2^m/(2^m + 1)
  return(beta)
}

# x: methylation data in beta values (columns:CpG sites; rows:people)
xtrain <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/norm_mvals.rds")%>%
  m2beta%>%t
xtest <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/wave3-final/w3.mvals.rds")%>%
  m2beta%>%t
inters_sites <- intersect(colnames(xtrain),colnames(xtest))
inters_450k <- intersect(inters_sites,probes)
methyl_train <- xtrain[,inters_450k]
methyl_test <- xtest[,inters_450k]

saveRDS(methyl_train,file="/Cluster_Filespace/Marioni_Group/Yufei/data/methyl_train.rds")
saveRDS(methyl_test,file="/Cluster_Filespace/Marioni_Group/Yufei/data/methyl_test.rds")

