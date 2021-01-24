library(psych)
int <- function(x){
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
}


setwd("D:/FangCloudV2/Zheng lab/Zheng lab共享资料/刘成林/ACS/ACS_code/explain/scfa_mol")

### all group explain
meta<-read.delim("scfa_mol_ratio.txt", row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
meta<-as.data.frame(meta)


### top75 pathway
pathway<-read.delim("pathway_rand.txt", row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
pathway_num<-ncol(pathway)

pathway<-as.data.frame(pathway)
pathway$ID<-rownames(pathway)
meta$ID<-rownames(meta)
### prevalence  >10 species
genus<-read.delim("species_rand.txt", row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
genus_num<-ncol(genus)
genus<-as.data.frame(genus)
genus_o<-genus

#####

genus$ID<-rownames(genus)
#####
genus<- merge(genus, pathway, by = "ID",all=FALSE)
##############
library(glmnet)

data_meta <- merge(meta, genus, by = "ID",all=FALSE)

daixie<-data_meta[,2:7]


bac<-data_meta[,8:ncol(data_meta)]


bac_notrans<-bac
for (i in 1:ncol(bac)){
  bac[,i]<-int(bac[,i])
}

for (i in 1:ncol(daixie)){
  daixie[,i]<-int(daixie[,i])
}


library(broom)
explatin_table<-matrix(NA,nrow = 100,ncol = ncol(daixie))

colnames(explatin_table)<-c(colnames(daixie))
for (j in 1:100){
  for (i in 1:ncol(daixie)){
    y<-as.matrix(daixie[,i])
    bac_log<-as.matrix(bac)
    cv_fit <- cv.glmnet(x=bac_log, y=y, nlambda = 1000,alpha = 1,family="gaussian")
    model_lasso_min <- glmnet(x=bac_log, y=y, alpha = 1, lambda=cv_fit$lambda.min,family="gaussian")
    explatin_table[j,i]<-model_lasso_min$dev.ratio
  }
}

### 100 result
write.csv(explatin_table,"all_explin_table_scfa_mol_100_af.csv",row.names = F)
### 100 mean
write.csv(t(as.data.frame(colMeans(explatin_table))),"all_explin_table_scfa_100_mol_mean_af.csv",row.names = F)
###############################################

####  case explain
meta<-read.delim("scfa_mol_ratio.txt", row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
meta<-as.data.frame(meta)
Case_index<-grep("A",rownames(meta))
meta<-meta[Case_index,]
### top75 pathway
pathway<-read.delim("pathway_rand.txt", row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
pathway_num<-ncol(pathway)

pathway<-as.data.frame(pathway)
pathway$ID<-rownames(pathway)
meta$ID<-rownames(meta)
### prevalence  >10 species
genus<-read.delim("species_rand.txt", row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
genus_num<-ncol(genus)
genus<-as.data.frame(genus)
genus_o<-genus

#####

genus$ID<-rownames(genus)
#####
genus<- merge(genus, pathway, by = "ID",all=FALSE)
##############
library(glmnet)

data_meta <- merge(meta, genus, by = "ID",all=FALSE)

daixie<-data_meta[,2:7]


bac<-data_meta[,8:ncol(data_meta)]


bac_notrans<-bac
for (i in 1:ncol(bac)){
  bac[,i]<-int(bac[,i])
}

for (i in 1:ncol(daixie)){
  daixie[,i]<-int(daixie[,i])
}


library(broom)
explatin_table<-matrix(NA,nrow = 100,ncol = ncol(daixie))

colnames(explatin_table)<-c(colnames(daixie))
for (j in 1:100){
  for (i in 1:ncol(daixie)){
    y<-as.matrix(daixie[,i])
    bac_log<-as.matrix(bac)
    cv_fit <- cv.glmnet(x=bac_log, y=y, nlambda = 1000,alpha = 1,family="gaussian")
    model_lasso_min <- glmnet(x=bac_log, y=y, alpha = 1, lambda=cv_fit$lambda.min,family="gaussian")
    explatin_table[j,i]<-model_lasso_min$dev.ratio
  }
}

### 100 result
write.csv(explatin_table,"Case_explin_table_scfa_mol_100_af.csv",row.names = F)
### 100 mean
write.csv(t(as.data.frame(colMeans(explatin_table))),"Case_explin_table_scfa_100_mol_mean_af.csv",row.names = F)


##################################### control

#### control explain
meta<-read.delim("scfa_mol_ratio.txt", row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
meta<-as.data.frame(meta)
Case_index<-grep("B",rownames(meta))
meta<-meta[Case_index,]
### top75 pathway
pathway<-read.delim("pathway_rand.txt", row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
pathway_num<-ncol(pathway)

pathway<-as.data.frame(pathway)
pathway$ID<-rownames(pathway)
meta$ID<-rownames(meta)
### prevalence  >10 species
genus<-read.delim("species_rand.txt", row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
genus_num<-ncol(genus)
genus<-as.data.frame(genus)
genus_o<-genus

#####

genus$ID<-rownames(genus)
#####
genus<- merge(genus, pathway, by = "ID",all=FALSE)
##############
library(glmnet)

data_meta <- merge(meta, genus, by = "ID",all=FALSE)

daixie<-data_meta[,2:7]


bac<-data_meta[,8:ncol(data_meta)]


bac_notrans<-bac
for (i in 1:ncol(bac)){
  bac[,i]<-int(bac[,i])
}

for (i in 1:ncol(daixie)){
  daixie[,i]<-int(daixie[,i])
}


library(broom)
explatin_table<-matrix(NA,nrow = 100,ncol = ncol(daixie))

colnames(explatin_table)<-c(colnames(daixie))
for (j in 1:100){
  for (i in 1:ncol(daixie)){
    y<-as.matrix(daixie[,i])
    bac_log<-as.matrix(bac)
    cv_fit <- cv.glmnet(x=bac_log, y=y, nlambda = 1000,alpha = 1,family="gaussian")
    model_lasso_min <- glmnet(x=bac_log, y=y, alpha = 1, lambda=cv_fit$lambda.min,family="gaussian")
    explatin_table[j,i]<-model_lasso_min$dev.ratio
  }
}

### 100 result
write.csv(explatin_table,"Control_explin_table_scfa_mol_100_af.csv",row.names = F)
### 100 mean
write.csv(t(as.data.frame(colMeans(explatin_table))),"Control_explin_table_scfa_100_mol_mean_af.csv",row.names = F)
