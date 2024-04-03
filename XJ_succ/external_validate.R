### import library & models
library(dplyr)
library(pROC)
library(ggplot2)
library(readxl)
source('/Users/pingyi/Desktop/XJ_succ/output/XJ_succ/used_function.R')
load(file = '/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/tssCA.Rdata')
load(file = '/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/tssWC.Rdata')
load(file = '/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/tssWIH.Rdata')

### import the external data set other CA types
external_clinic <- read_excel("~/Desktop/XJ/XJ_paper/data/clinic/external_clinic.xlsx")
matrix_tss_external <- read.csv("~/Desktop/XJ/XJ_paper/data/rawdata/matrix_tss_external.csv",row.names = 1)
dim(matrix_tss_external) # 93701   246
### import the data for building model
data_model <- read_excel("~/Desktop/XJ_succ/output/APaper/sup1.xlsx",row.names = 1)
data_model <- data_model[,2:ncol(data_model)]
colnames(data_model) <- data_model[1,]
data_model <- data_model[2:nrow(data_model),]
### import the external data set of W & C
pred_ex_all <- read.csv("~/Desktop/pred_ex_all.csv", header=FALSE)
pred_ex_all$patient <- apply(data.frame(pred_ex_all$V1),1,
                             function(x){return(strsplit(x,'[.]')[[1]][1])})
### pred_ex_all samples not used for building model
pred_ex_all$inmodel_set <- apply(data.frame(pred_ex_all$patient),1,
                                 function(x){return(x %in% data_model$patient )})
pred_ex_all <- pred_ex_all[which(pred_ex_all$inmodel_set  == 'FALSE'),]
dim(pred_ex_all) #  203  12
### rm duplicated patients
pred_ex_all$dupornot <- duplicated(pred_ex_all$patient)
pred_ex_all <- pred_ex_all[which(pred_ex_all$dupornot  == 'FALSE'),]
dim(pred_ex_all) #  176  12
### rm NA rows
pred_ex_all <- pred_ex_all[!is.na(pred_ex_all$V2),]
dim(pred_ex_all) #  175  13
### plot the ROC curve -- calculate the auc score
combinde_auc <- roc(pred_ex_all$V7,pred_ex_all$V2,type="prob") # 0.9071
### import the tss matrix of this validation cohort
tssmatrix_nova <- read.delim("~/Desktop/XJ_new/datav2/tss/tssmatrix_nova.txt")
tssmatrix_T7 <- read.delim("~/Desktop/XJ_new/datav2/tss/tssmatrix_T7.txt")

extssmatrix <- tssmatrix_nova[,colnames(tssmatrix_nova) %in% pred_ex_all$V1]
dim(extssmatrix) # 93701   129
extssmatrix2 <- matrix_tss_external[,colnames(matrix_tss_external) %in% pred_ex_all$V1]
dim(extssmatrix2) # 93701   140
extssmatrix3 <- cbind(extssmatrix,extssmatrix2)
unique(colnames(extssmatrix3))

extssmatrix3 <- extssmatrix3[,!duplicated(colnames(extssmatrix3))]
dim(extssmatrix3) #  93701   174
#setdiff(pred_ex_all$V1,colnames(extssmatrix3)) # "WITT0P0015.1.TWN1"
### site number 
site_used <- rownames(data.frame(modeltss_CAnonCA$finalModel$importance))
matrix_tss_external <- extssmatrix3[site_used,]
### z-score treat as the training_range_site
matrix_tss_externalt <- data.frame(t(matrix_tss_external))
external_zscore <- newmatrix_z_norm(new_matrix=matrix_tss_externalt,training_range_site)
dim(external_zscore) #  174 90729
### choose features as train_topbins did
dim(train_topbins) #  338 301
site_trainN <-300
external_zscore_feat <- external_zscore[,colnames(train_topbins)[1:site_trainN]]
dim(external_zscore_feat) #  174 300
### pred with the model 'modeltss_CAnonCA'
pred_ex <- predict(modeltss_CAnonCA,newdata = external_zscore_feat,type = "prob")
pred_ex$sample = rownames(pred_ex)
####################################################################################
####################### check the result of 'modeltss_WC' ##########################
####################################################################################
### site number 
site_used_WC <- rownames(data.frame(modeltss_WC$finalModel$importance))
matrix_tss_externalWC <- external_zscore[,site_used_WC]
dim(matrix_tss_externalWC) # 174 300
### z-score treat as the training_range_site
external_zscoreWC <- matrix_tss_externalWC
dim(external_zscoreWC) #   174 300
### pred with the model 'modeltss_CAnonCA'
pred_exWC <- predict(modeltss_WC,newdata = external_zscoreWC,type = "prob")
pred_exWC$sample = rownames(pred_exWC)
pred_exWC$source = apply(data.frame(pred_exWC$sample),1,
                         function(x){ return(strsplit(x,'')[[1]][1])})

pred_exWC <- pred_exWC[which(pred_exWC$source %in% c('W','C')),] # 142
pred_exWC$true = 1
pred_exWC[which(pred_exWC$source == 'C'),'true'] = 0

### plot the ROC curve -- calculate the auc score
combinde_auc_WC <- roc(pred_exWC$true,pred_exWC$W,type="prob") # 0.9426
