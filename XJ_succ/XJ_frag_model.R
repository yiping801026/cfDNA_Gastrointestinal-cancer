### This script is used to build the models based on fragmentation profile 
### Three models have the goals to 1. Early cancer screening model (CA vs non-CA) 
### 2. Early cancer screening model (Healthy vs Gastritis)  
### 3. Tissue-of-origin localization model (Gastric cancer vs  Colorectal cancer)
### use the pre-treated norm dataset of fragmentation profile 'df_ratio_norm_sample.csv'
### 'df_ratio_norm_sample.csv' is treated by the script ''

### import library
library(readxl)
library(dplyr)
library(factoextra) 
library(ggplot2)
library(pROC)
library(stats)
library(caret)

### import inhouse functions
source('/Users/pingyi/Desktop/XJ_succ/output/XJ_succ/used_function.R')

####################################################################################
############################        CA nonCA              ##########################
####################################################################################

### import the dataset
### import the clinic information 
model_clinic <- read_excel("/Users/pingyi/Desktop/XJ_succ/output/data/model_clinic2.xlsx")
matrix_frag <- read.csv("/Users/pingyi/Desktop/XJ_succ/output/data/df_ratio_norm_sample.csv", row.names=1)
#double_Healthy_clinic <- model_clinic
matrix_frag_model=matrix_frag[,model_clinic$sample]
dim(matrix_tss_model)
### pretreat 
#matrix_tss_model <- matrix_tss_model[,-1]
matrix_tss_model_pred <- pre_treat_rm0pNA(matrix_tss_model)
#matrix_tss_model_nzv <- nearZeroVar_for_matrix(matrix_tss_model_pred)

dim(matrix_tss_model_pred) #   451 90717
ll ='none'
### change >10 into 10
if(ll == 'none'){
  matrix_tss_model_nzv <- matrix_tss_model_pred
}else  {
  print(c('bigger tss scores will be change into',as.numeric(args[2])))
  cut_big = as.numeric(args[2])
  matrix_tss_model_nzv[matrix_tss_model_nzv>as.numeric(cut_big)] = as.numeric(cut_big)
  
}
#matrix_tss_model_nzv[matrix_tss_model_nzv>10] = 10

### 75% for training 25% for testing
set.seed(123)
inTrain=createDataPartition(y=model_clinic$Canon,p=0.75,list=FALSE)
training=matrix_tss_model_nzv[inTrain,]
testing=matrix_tss_model_nzv[-inTrain,]
### save the group of training && testing 
trainLabels <- model_clinic[inTrain,4]
testLabels <- model_clinic[-inTrain,4]
test_clinic_matrix<- model_clinic[-inTrain,]
dim(training) #  339 90717
dim(testing)  # 112 90717


training_range_sample <- range_cal(matrix = training,roworcol = 1)
range(training_range_sample$median) # 0.00000 1.17108  

training_range_site <- range_cal(matrix = training,roworcol = 2)
range(training_range_site$median) # 0.0000 117.5705
training_zscore <- training
testing_zscore <- testing

alltss_zscore <- rbind(training_zscore,testing_zscore)
################### rm batch
DH_zscore <- alltss_zscore[double_Healthy_clinic$sample,]
DH_zscore_t = data.frame(t(DH_zscore))

double_Healthy_clinic$batch <- 'one'
double_Healthy_clinic[which(double_Healthy_clinic$platform == 'T7'),ncol(double_Healthy_clinic)] = 'two'

tss_z_dhpbatch <- p_cal(trainTransformed = DH_zscore,
                        trainLabels =  double_Healthy_clinic$batch,
                        trainLabels_sig = 'one',
                        paired = F) #    


### remove the p<0.05 sites to choose features 

tss_z_dhpbatchcut <- tss_z_dhpbatch[which(tss_z_dhpbatch[,1] > as.numeric(0.01) ),]
dim(tss_z_dhpbatchcut) #   470   2
DH_zscore_rmbatch <- DH_zscore[,tss_z_dhpbatchcut$V2]
dim(DH_zscore_rmbatch) #    450 43352
DH_zscore_rmbatch_t = data.frame(t(DH_zscore_rmbatch))

training_zscore_rmbatch <- training_zscore[,tss_z_dhpbatchcut$V2]
dim(training_zscore_rmbatch) #   339 30547
training_zscore_rmbatch_t = data.frame(t(training_zscore_rmbatch))
range(training_zscore_rmbatch_t) # -4.44812 4817.93667
frag_norm_rmbatch<- DH_zscore_rmbatch

platform_train = model_clinic[inTrain,'platform']


group_train = model_clinic[inTrain,'Canon']

tss_z_pCA <- p_cal(trainTransformed = training_zscore_rmbatch,
                   trainLabels = group_train$Canon,
                   trainLabels_sig = 'CA',paired = FALSE)# 2.30879442184768e-10   0.999946126323105  

dim(tss_z_pCA) # 43352     2
tss_z_pCAcut <- tss_z_pCA[which(tss_z_pCA[,1] < as.numeric(1) ),]
dim(tss_z_pCAcut)  # 193   2
### choose top 300
tss_z_pCAcut <- tss_z_pCAcut[order(tss_z_pCAcut$tout),]
tss_z_pCAcut <- tss_z_pCAcut[1:as.numeric(300),]

###
dim(tss_z_pCAcut) #  3239    2
frag_CAnon_grouppsites <- tss_z_pCAcut
### build the model based on features:tss_z_pCAcut$V2

train_topbins<-as.data.frame(training_zscore_rmbatch)[,tss_z_pCAcut$V2]
dim(train_topbins) #  339 3239
test_topbins<-as.data.frame(testing_zscore)[,tss_z_pCAcut$V2]
dim(test_topbins) # 112 3239

train_topbins$group <- factor(trainLabels$Canon,levels = c('CA','nonCA'))
test_topbins$group <- factor(testLabels$Canon,levels = c('CA','nonCA'))

fragCAnontrain_topbins <- train_topbins
fragCAnontest_topbins <- test_topbins





### build model based on rf
set.seed(123) 
ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 3,
                    allowParallel = TRUE,classProbs = TRUE,savePredictions = TRUE,search='grid',summaryFunction = twoClassSummary)
set.seed(123) 
ctrl$sampling <- "up"
tunegrid <- expand.grid(mtry = c(1:5))
set.seed(123)
### 0.9027833   0.9027
rf_default1 <- train(group~., 
                     data=train_topbins, 
                     method='rf', 
                     #metric='Accuracy',  #Metric compare model is Accuracy
                     metric = "ROC",
                     tuneGrid=tunegrid,
                     ntree=c(500),
                     trControl=ctrl)

modelfrag_CAnonCA <- rf_default1

rf_default1
pred <- predict(rf_default1,newdata = test_topbins,type = "prob")
pred$sample = rownames(pred)
pred <- left_join(pred,model_clinic,by='sample')
pred$res = 1
pred[which(pred$CA<0.5),'res'] = 0 

pred$true = 1
pred[which(pred$Canon != 'CA'),'true'] = 0 

pred$equ = 1
pred[which(pred$res != pred$true),'equ'] = 0 
View(pred)

frag_CAnopred <- pred
test_obs <- data.frame(test_topbins$group)
test_obs$obs_prob <- 'm'
test_obs[which(test_obs$test_topbins.group=='nonCA'),]$obs_prob = 0
test_obs[which(test_obs$test_topbins.group=='CA'),]$obs_prob = 1
roc_data <- roc(test_obs$obs_prob, pred[,1],type="prob")
print(roc_data)

frag_CAnoroc_data <- roc_data
save(frag_norm_rmbatch,frag_CAnon_grouppsites,fragCAnontrain_topbins,
     fragCAnontest_topbins,modelfrag_CAnonCA,frag_CAnopred,frag_CAnoroc_data,
     file='/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/fragCA.Rdata')


####################################################################################
############################         H vs  inflammation   ##########################
####################################################################################
load('/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/fragCA.Rdata')
### use the data rm batch (DH_zscore_rmbatch_t) clinic (model_clinic)
tssrmbatch <- frag_norm_rmbatch
model_clinic <- read_excel("~/Desktop/XJ/XJ_paper/data/clinic/model_clinic2.xlsx")
###only WI and H included
model_clinic_WIH <- model_clinic[which(model_clinic$disease !='T'),]
tssrmbatch <- tssrmbatch[model_clinic_WIH$sample,]

dim(tssrmbatch)
### 75% for training 25% for testing
set.seed(123)
inTrain=createDataPartition(y=model_clinic_WIH$disease,p=0.75,list=FALSE)
training=tssrmbatch[inTrain,]
testing=tssrmbatch[-inTrain,]
### save the group of training && testing 
trainLabels <- model_clinic_WIH[inTrain,'disease']
testLabels <- model_clinic_WIH[-inTrain,'disease']
test_clinic_matrix<- model_clinic_WIH[-inTrain,]
dim(training) # 197 43352
dim(testing)  #  65 43352

### calculate the p-value between groups

group_train = model_clinic[inTrain,'disease']

tss_z_pCA <- p_cal(trainTransformed = training,
                   trainLabels = group_train$disease,
                   trainLabels_sig = 'I',paired = FALSE)# 2.30879442184768e-10   0.999946126323105  

dim(tss_z_pCA) #
tss_z_pCAcut <- tss_z_pCA[which(tss_z_pCA[,1] < as.numeric(1) ),]
dim(tss_z_pCAcut)  #3262    2
### choose top 300
tss_z_pCAcut <- tss_z_pCAcut[order(tss_z_pCAcut$tout),]
tss_z_pCAcut <- tss_z_pCAcut[1:as.numeric(300),]

###
dim(tss_z_pCAcut) #  300   2
### build the model based on features:tss_z_pCAcut$V2

train_topbins<-as.data.frame(training)[,tss_z_pCAcut$V2]
dim(train_topbins) #  339 3239
test_topbins<-as.data.frame(testing)[,tss_z_pCAcut$V2]
dim(test_topbins) # 112 3239

train_topbins$group <- factor(trainLabels$disease,levels = c('I','H'))
test_topbins$group <- factor(testLabels$disease,levels = c('I','H'))


### build model based on rf
set.seed(123) 
ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 3,
                    allowParallel = TRUE,classProbs = TRUE,savePredictions = TRUE,search='grid',summaryFunction = twoClassSummary)
set.seed(123) 
ctrl$sampling <- "up"
tunegrid <- expand.grid(mtry = c(1:3))
set.seed(123)
### 0.9027833   0.9027
modelfrag_WIH <- train(group~., 
                      data=train_topbins, 
                      method='rf', 
                      #metric='Accuracy',  #Metric compare model is Accuracy
                      metric = "ROC",
                      tuneGrid=tunegrid,
                      ntree=c(500),
                      trControl=ctrl)

rf_default1 <- modelfrag_WIH
rf_default1
#tss_CAnonCA1000 <- rf_default1
pred <- predict(rf_default1,newdata = test_topbins,type = "prob")
pred$sample = rownames(pred)
pred <- left_join(pred,model_clinic,by='sample')
pred$res = 1
pred[which(pred$I<0.5),'res'] = 0 

pred$true = 1
pred[which(pred$disease != 'I'),'true'] = 0 

pred$equ = 1
pred[which(pred$res != pred$true),'equ'] = 0 
View(pred)


test_obs <- data.frame(test_topbins$group)
test_obs$obs_prob <- 'm'
test_obs[which(test_obs$test_topbins.group=='H'),]$obs_prob = 0
test_obs[which(test_obs$test_topbins.group=='I'),]$obs_prob = 1
roc_data <- roc(test_obs$obs_prob, pred[,1],type="prob")
print(roc_data)
#print(c(args[1],args[2],args[3],args[4]))


pred_fragWIH <- pred
roc_data_fragWIH <- roc_data
train_topbins_fragWIH <- train_topbins 
test_topbins_fragWIH<- test_topbins


save(train_topbins_fragWIH,test_topbins_fragWIH,
     modelfrag_WIH,pred_fragWIH,roc_data_fragWIH,
     file='/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/fragWIH.Rdata')


####################################################################################
############################    localization: G vs  C ##############################
####################################################################################
load('/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/fragCA.Rdata')
### use the data rm batch (DH_zscore_rmbatch_t) clinic (model_clinic)
tssrmbatch <- frag_norm_rmbatch
model_clinic <- read_excel("~/Desktop/XJ/XJ_paper/data/clinic/model_clinic2.xlsx")
###only WI and H included
model_clinic_WC <- model_clinic[which(model_clinic$disease =='T'),]
tssrmbatch <- tssrmbatch[model_clinic_WC$sample,]

dim(tssrmbatch) #188 43352
### 75% for training 25% for testing
set.seed(123)
inTrain=createDataPartition(y=model_clinic_WC$disease,p=0.75,list=FALSE)
training=tssrmbatch[inTrain,]
testing=tssrmbatch[-inTrain,]
### save the group of training && testing 
trainLabels <- model_clinic_WC[inTrain,'source']
testLabels <- model_clinic_WC[-inTrain,'source']
test_clinic_matrix<- model_clinic_WC[-inTrain,]
dim(training) #  141 43352
dim(testing)  #  65 43352

### calculate the p-value between groups

group_train = model_clinic_WC[inTrain,'source']

tss_z_pCA <- p_cal(trainTransformed = training,
                   trainLabels = group_train$source,
                   trainLabels_sig = 'W',paired = FALSE)# 2.30879442184768e-10   0.999946126323105  

dim(tss_z_pCA) #43352     2
tss_z_pCAcut <- tss_z_pCA[which(tss_z_pCA[,1] < as.numeric(0.05) ),]
dim(tss_z_pCAcut)  # 1201    2
### choose top 300
tss_z_pCAcut <- tss_z_pCAcut[order(tss_z_pCAcut$tout),]
#tss_z_pCAcut <- tss_z_pCAcut[1:as.numeric(300),]

###
dim(tss_z_pCAcut) #  300   2
### build the model based on features:tss_z_pCAcut$V2

train_topbins<-as.data.frame(training)[,tss_z_pCAcut$V2]
dim(train_topbins) #  339 3239
test_topbins<-as.data.frame(testing)[,tss_z_pCAcut$V2]
dim(test_topbins) # 112 3239

train_topbins$group <- factor(trainLabels$source,levels = c('W','C'))
test_topbins$group <- factor(testLabels$source,levels = c('W','C'))


### build model based on rf
set.seed(123) 
ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 3,
                    allowParallel = TRUE,classProbs = TRUE,savePredictions = TRUE,search='grid',summaryFunction = twoClassSummary)
set.seed(123) 
ctrl$sampling <- "up"
tunegrid <- expand.grid(mtry = c(1:5))
set.seed(123)
### 0.9027833   0.9027
modelfrag_WC <- train(group~., 
                     data=train_topbins, 
                     method='rf', 
                     #metric='Accuracy',  #Metric compare model is Accuracy
                     metric = "ROC",
                     tuneGrid=tunegrid,
                     ntree=c(500),
                     trControl=ctrl)

rf_default1 <- modelfrag_WC
rf_default1
#tss_CAnonCA1000 <- rf_default1
pred <- predict(rf_default1,newdata = test_topbins,type = "prob")
pred$sample = rownames(pred)
pred <- left_join(pred,model_clinic,by='sample')
pred$res = 1
pred[which(pred$W<0.5),'res'] = 0 

pred$true = 1
pred[which(pred$source != 'W'),'true'] = 0 

pred$equ = 1
pred[which(pred$res != pred$true),'equ'] = 0 
View(pred)


test_obs <- data.frame(test_topbins$group)
test_obs$obs_prob <- 'm'
test_obs[which(test_obs$test_topbins.group=='C'),]$obs_prob = 0
test_obs[which(test_obs$test_topbins.group=='W'),]$obs_prob = 1
roc_data <- roc(test_obs$obs_prob, pred[,1],type="prob")
print(roc_data)
#print(c(args[1],args[2],args[3],args[4]))


pred_fragWC <- pred
roc_data_fragWC <- roc_data
train_topbins_fragWC<- train_topbins 
test_topbins_fragWC<- test_topbins


save(train_topbins_fragWC,test_topbins_fragWC,
     modelfrag_WC,pred_fragWC,roc_data_fragWC,
     file='/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/fragWC.Rdata')
/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/fragWIH.Rdata












