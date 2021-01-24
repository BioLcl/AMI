
int <- function(x){
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
}
################################## all group

setwd("D:/FangCloudV2/Zheng lab/Zheng lab共享资料/刘成林/ACS/ACS_code/explain/TMAO")


meta<-read.delim("TMAO_select.txt", row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
meta_num<-ncol(meta)


meta<-as.data.frame(meta)
# all_index<-grep("A",rownames(meta))
# # length(all_index)
# 
# meta<-meta[all_index,]


pathway<-read.delim("pathway_rand.txt", row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
pathway_num<-ncol(pathway)

pathway<-as.data.frame(pathway)
pathway$ID<-rownames(pathway)
meta$ID<-rownames(meta)

genus<-read.delim("species_rand.txt", row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
genus_num<-ncol(genus)


genus<-as.data.frame(genus)





#####
genus_o<-genus

#####

genus$ID<-rownames(genus)
#####
genus<- merge(genus, pathway, by = "ID",all=FALSE)
##############
library(glmnet)

data_meta <- merge(meta, genus, by = "ID",all=FALSE)
###

###
daixie<-data_meta[,2:(meta_num+1)]


bac<-data_meta[,(meta_num+2):ncol(data_meta)]
bac<-bac[,colSums(bac)!=0]

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
    y_temp<-y
    colnames(y_temp)<-colnames(daixie)[i]
    y<-y[!is.na(y_temp)]
    bac_log<-as.matrix(bac)[!is.na(y_temp),]
    
    cv_fit <- cv.glmnet(x=bac_log, y=y, nlambda = 1000,alpha = 1,family="gaussian")
    model_lasso_min <- glmnet(x=bac_log, y=y, alpha = 1, lambda=cv_fit$lambda.min,family="gaussian")
    explatin_table[j,i]<-model_lasso_min$dev.ratio
  }
}


write.csv(explatin_table,"all_TMAO_explin_table_100_af.csv",row.names = F)
write.csv(t(as.data.frame(colMeans(explatin_table))),"all_TMAO_explin_table_100_mean_af.csv",row.names = F)


############## case group
meta<-read.delim("TMAO_select.txt", row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
meta_num<-ncol(meta)


meta<-as.data.frame(meta)
all_index<-grep("A",rownames(meta))

meta<-meta[all_index,]


pathway<-read.delim("pathway_rand.txt", row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
pathway_num<-ncol(pathway)

pathway<-as.data.frame(pathway)
pathway$ID<-rownames(pathway)
meta$ID<-rownames(meta)

genus<-read.delim("species_rand.txt", row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
genus_num<-ncol(genus)


genus<-as.data.frame(genus)





#####
genus_o<-genus

#####

genus$ID<-rownames(genus)
#####
genus<- merge(genus, pathway, by = "ID",all=FALSE)
##############
library(glmnet)

data_meta <- merge(meta, genus, by = "ID",all=FALSE)
###

###
daixie<-data_meta[,2:(meta_num+1)]


bac<-data_meta[,(meta_num+2):ncol(data_meta)]
bac<-bac[,colSums(bac)!=0]

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
    y_temp<-y
    colnames(y_temp)<-colnames(daixie)[i]
    y<-y[!is.na(y_temp)]
    bac_log<-as.matrix(bac)[!is.na(y_temp),]
    
    cv_fit <- cv.glmnet(x=bac_log, y=y, nlambda = 1000,alpha = 1,family="gaussian")
    model_lasso_min <- glmnet(x=bac_log, y=y, alpha = 1, lambda=cv_fit$lambda.min,family="gaussian")
    explatin_table[j,i]<-model_lasso_min$dev.ratio
  }
}


write.csv(explatin_table,"case_TMAO_explin_table_100_af.csv",row.names = F)
write.csv(t(as.data.frame(colMeans(explatin_table))),"case_TMAO_explin_table_100_mean_af.csv",row.names = F)


############## control group
meta<-read.delim("TMAO_select.txt", row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
meta_num<-ncol(meta)


meta<-as.data.frame(meta)
all_index<-grep("B",rownames(meta))

meta<-meta[all_index,]


pathway<-read.delim("pathway_rand.txt", row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
pathway_num<-ncol(pathway)

pathway<-as.data.frame(pathway)
pathway$ID<-rownames(pathway)
meta$ID<-rownames(meta)

genus<-read.delim("species_rand.txt", row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
genus_num<-ncol(genus)


genus<-as.data.frame(genus)





#####
genus_o<-genus

#####

genus$ID<-rownames(genus)
#####
genus<- merge(genus, pathway, by = "ID",all=FALSE)
##############
library(glmnet)

data_meta <- merge(meta, genus, by = "ID",all=FALSE)
###

###
daixie<-data_meta[,2:(meta_num+1)]


bac<-data_meta[,(meta_num+2):ncol(data_meta)]
bac<-bac[,colSums(bac)!=0]

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
    y_temp<-y
    colnames(y_temp)<-colnames(daixie)[i]
    y<-y[!is.na(y_temp)]
    bac_log<-as.matrix(bac)[!is.na(y_temp),]
    
    cv_fit <- cv.glmnet(x=bac_log, y=y, nlambda = 1000,alpha = 1,family="gaussian")
    model_lasso_min <- glmnet(x=bac_log, y=y, alpha = 1, lambda=cv_fit$lambda.min,family="gaussian")
    explatin_table[j,i]<-model_lasso_min$dev.ratio
  }
}


write.csv(explatin_table,"control_TMAO_explin_table_100_af.csv",row.names = F)
write.csv(t(as.data.frame(colMeans(explatin_table))),"control_TMAO_explin_table_100_mean_af.csv",row.names = F)





