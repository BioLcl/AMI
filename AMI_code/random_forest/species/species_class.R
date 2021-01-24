library("tidyverse")
library("caret")
library("pROC")
library("randomForest")
library("data.table")
library("rfUtilities")
library("plyr")
library("e1071")
metadetail<-1
split_number=100
bac<-read.delim("species_rand.txt",header = T,stringsAsFactors = F,check.names = F,row.names = 1)
# daixie<-read.delim("groupdata.txt",header = T,stringsAsFactors = F,check.names = F,row.names = 1)
# daixie$class<-"Case"
# daixie[113:nrow(daixie),1]<-"Control"
# write.csv(daixie,"groupdata.csv")
daixie<-read.csv("groupdata.csv",header = T,stringsAsFactors = F,check.names = F,row.names = 1)


#bac<-as.data.frame(t(bac))


bac<-bac[match(rownames(daixie),rownames(bac)),]


remove_rare <- function( table , cutoff_pro ) {
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * nrow(table) )
  for ( i in 1:ncol(table) ) {
    row_nonzero <- length( which( table[  , i]  > 0 ) )
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  return( table [ , row2keep, drop=F ])
}

#bac<-remove_rare(bac,0.1)
#############
bac$names<-rownames(bac)
daixie$names<-rownames(daixie)
dataall<-merge(daixie,bac,by="names",all=F)
dataall<-dataall[,-1]
dataset<-dataall[,c(1,2:ncol(dataall))]
daixieID<-colnames(dataset)[1]
dir<-paste0("./",daixieID)
dir.create(dir)
colnames(dataset)[1]<-"classes"

#preProcValues <- preProcess(dataset, method = "range")
#dataset <- predict(preProcValues, dataset)

#pipeline <- function(dataset,daixieID){
  # Create vectors to save cv and test AUC values for every data-split
  results_total <-  data.frame()
  test_aucs <- c()
  cv_aucs <- c()
  dataTransformed<-dataset
  # We are doing the pre-processing to the full dataset and then splitting 80-20
  # Scale all features between 0-1
  # preProcValues <- preProcess(dataset, method = "range")
  # dataTransformed <- predict(preProcValues, dataset)
  # 
  # Do the 80-20 data-split
  # Stratified data partitioning %80 training - %20 testing
  inTraining <- createDataPartition(dataTransformed$classes, p = .70, list = FALSE)
  trainTransformed <- dataTransformed[ inTraining,]
  testTransformed  <- dataTransformed[-inTraining,]

preProcValues <- preProcess(trainTransformed, method = "range")
trainTransformed<-predict(preProcValues,trainTransformed)
testTransformed<-predict(preProcValues,testTransformed)
  # remove columns that only appear within one or fewer samples of the training set. These are
  # likely to be all zero and will not enter into the model
  # frequent <- names(which(apply(dataTransformed[, -1] > 0, 2, sum) > 1))
  # trainTransformed <- trainTransformed %>% select(classes, frequent)
  # testTransformed <- testTransformed %>% select(classes, frequent)
  

#    rfeControls_rf <- rfeControl(
#    functions = rfFuncs,
#    method = 'cv',
#    repeats = 5)
#
#    trainTransformed[,1]<-as.factor(trainTransformed[,1])
#  fs_nb <- rfe(x = trainTransformed[,-1],
#               y = trainTransformed[,1],
#               sizes = seq(1,ncol(dataTransformed)-1,3),
#               rfeControl = rfeControls_rf)
#  fs_nb$optVariables
#  write.csv(fs_nb$optVariables,paste("./",daixieID,"/","select_feature.csv",sep=""))
#  selectbac<-bac[,fs_nb$optVariables]
#  write.csv(selectbac,paste("./",daixieID,"/","select_bac.csv",sep=""))
#3  testTransformed<-testTransformed[,c("class",fs_nb$optVariables)]
#  trainTransformed<-trainTransformed[,c("class",fs_nb$optVariables)]

  write.csv(trainTransformed,paste("./",daixieID,"/","train_data.csv",sep=""))
  write.csv(testTransformed,paste("./",daixieID,"/","test_data.csv",sep=""))

 
  
  
  #######################
  n_features <- ncol(trainTransformed) - 1
  if(n_features > 20000) n_features <- 20000
  
  if(n_features < 19){ mtry <- 1:6
  } else { mtry <- floor(seq(1, n_features/3, length=6)) }
  
  mtry <- mtry[mtry <= n_features]
  
  # cv index to make sure the internal 5-folds are stratified for diagnosis classes and also resampled 100 times.
  # 100 repeat internally is necessary to get robust readings of hyperparameter setting performance
  folds <- 10
  cvIndex <- createMultiFolds(factor(trainTransformed$classes), folds, times=10) #returnTrain = T default for multifolds
  
  cv <- trainControl(method="repeatedcv",
                     number= folds,
                     index = cvIndex,
                     returnResamp="final",
                     classProbs=TRUE,
                     summaryFunction=twoClassSummary,
                     indexFinal=NULL,
                     savePredictions = TRUE)
  
  grid <-  expand.grid(mtry = mtry)
  
  # Train the model
  set.seed(2)
  trained_model <-  train(classes ~ .,
                          data=trainTransformed,
                          method = "rf",
                          importance=TRUE,
                          trControl = cv,
                          metric = "ROC",
                          tuneGrid = grid,
                          ntree=1000)
  imp= as.data.frame(trained_model$finalModel$importance)
  imp = imp[order(imp[,3],decreasing = T),]
  write.table(imp,paste("./",daixieID,"/","importance.csv",sep=""))
  # Mean AUC value over repeats of the best cost parameter during training
  cv_auc <- getTrainPerf(trained_model)$TrainROC
  
  # Predict on the test set and get predicted probabilities
  rpartProbs <- predict(trained_model, testTransformed, type="prob")
  test_roc <- roc(ifelse(testTransformed$classes == "Case", 1, 0),
                  rpartProbs[[2]])
  test_auc <- test_roc$auc
  pdf(file ="./test_auc.pdf",width =6,height =6)
  plot(test_roc, print.auc = TRUE, auc.polygon = TRUE, legacy.axes = TRUE, 
       grid = c(0.1, 0.2), grid.col = c("green", "red"), max.auc.polygon = TRUE,  
       auc.polygon.col = "skyblue", print.thres = TRUE, 
       main = "random forest") 
  dev.off()
  # Save all the test AUCs over iterations in test_aucs
  test_aucs <- c(test_aucs, test_auc)
  
  # Cross-validation mean AUC value
  # Save all the cv meanAUCs over iterations in cv_aucs
  cv_aucs <- c(cv_aucs, cv_auc)
  
  # Save all results of hyper-parameters and their corresponding meanAUCs for each iteration
  results_individual <- trained_model$results
  results_total <- rbind(results_total, results_individual)
  
  results <- list(cv_aucs, test_aucs, results_total)
#  return(results)
#}

#results <- pipeline(dataset,daixieID)
aucs <- matrix(c(results[[1]], results[[2]]), ncol=2)
# Convert to dataframe and add a column noting the model name
aucs_dataframe <- data.frame(aucs)
colnames(aucs_dataframe)<-c("cv_aucs","test_aucs")
write_csv(aucs_dataframe,path=paste0(dir, "/optimum_mtry.", split_number, ".csv"))
# ------------------------------------------------------------------
# Save all tunes from 100 data splits and corresponding AUCs
all_results <- results[3]
# Convert to dataframe and add a column noting the model name
dataframe <- data.frame(all_results) %>%
  write_csv(path=paste0(dir, "/all_mtry.", split_number, ".csv"))
auc_data<-cbind(test_roc$sensitivities,test_roc$specificities)
colnames(auc_data)<-c("sensitivities","specificities")
write.table(auc_data,"auc_species_data.txt",sep="\t",quote=F,row.names=F)
write.table(ci.auc(test_roc),auc_ci_data.txt")
save.image("class_spe.RData")

