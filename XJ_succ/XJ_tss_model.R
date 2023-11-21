### This script is used to build the model based on tss dataset
### There are three tasks for this program:
# 1. to distinguish the cancer versus non-cancer samples (CA vs nonCA), nonCA 
# samples include healthy individuals and samples with inflammation
# 2. to distinguish the healthy and gastritis individuals (H vs  inflammation)
# 3. to distinguish cancer localization (G vs C), gastric vs colorectal cancer

### import library
library(readxl)
library(dplyr)
library(factoextra) 
library(ggplot2)
library(pROC)
library(stats)
library(caret)

####################################################################################
############################        CA nonCA              ##########################
####################################################################################

### import inhouse functions
source('./used_function.R')

### import dataset
model_clinic <- read_excel("/Users/pingyi/Desktop/XJ_succ/output/data/model_clinic2.xlsx")
double_Healthy_clinic <- model_clinic
matrix_tss <- read.csv("/Users/pingyi/Desktop/XJ_succ/output/data/matrix_tss.csv", row.names=1)
dim(matrix_tss)

matrix_tss_model <- matrix_tss[,model_clinic$sample]
dim(matrix_tss_model)
### pretreat 
#matrix_tss_model <- matrix_tss_model[,-1]
matrix_tss_model_pred <- pre_treat_rm0pNA(matrix_tss_model)
#matrix_tss_model_nzv <- nearZeroVar_for_matrix(matrix_tss_model_pred)

dim(matrix_tss_model_pred) #   451 90717
ll ='none'
### 
if(ll == 'none'){
  matrix_tss_model_nzv <- matrix_tss_model_pred
}else  {
  print(c('bigger tss scores will be change into',as.numeric(10)))
  cut_big = as.numeric(10)
  matrix_tss_model_nzv[matrix_tss_model_nzv>as.numeric(cut_big)] = as.numeric(cut_big)
  
}

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

### do the normalize for the sites -- use the z-score && median & MAD and ignore the outliers
training_range_site <- training_range_site[which(training_range_site$mad != 0),]
site_used <- rownames(training_range_site[which(training_range_site$mad != 0),])

training_rmmad <- training[,site_used]

training_zscore <- data.frame(apply(data.frame(training_rmmad), 2, custom_z_score))
training_zscore_range <- range_cal(matrix = training_zscore,roworcol = 1)
range(training_zscore_range$median) # -1.4207804  0.8213466
dim(training_zscore) #  248 90808

testing_zscore <- newmatrix_z_norm(new_matrix=testing,training_range_site)
dim(testing_zscore) # 82 90808
testing_zscore_range <- range_cal(matrix = testing_zscore,roworcol = 1)
range(testing_zscore_range$median) # -0.6744908  0.5522655




alltss_zscore <- rbind(training_zscore,testing_zscore)
#double_Healthy_clinic
DH_zscore <- alltss_zscore[double_Healthy_clinic$sample,]
DH_zscore_t = data.frame(t(DH_zscore))



double_Healthy_clinic$batch <- 'one'
double_Healthy_clinic[which(double_Healthy_clinic$platform == 'T7'),ncol(double_Healthy_clinic)] = 'two'

tss_z_dhpbatch <- p_cal(trainTransformed = DH_zscore,
                        trainLabels =  double_Healthy_clinic$batch,
                        trainLabels_sig = 'one',
                        paired = F) #    


### remove the p<0.05 sites to choose features 

tss_z_dhpbatchcut <- tss_z_dhpbatch[which(tss_z_dhpbatch[,1] > as.numeric(0.05) ),]
dim(tss_z_dhpbatchcut) #   30547     2
DH_zscore_rmbatch <- DH_zscore[,tss_z_dhpbatchcut$V2]
dim(DH_zscore_rmbatch) #    450 43352
DH_zscore_rmbatch_t = data.frame(t(DH_zscore_rmbatch))

training_zscore_rmbatch <- training_zscore[,tss_z_dhpbatchcut$V2]
dim(training_zscore_rmbatch) #   339 30547
training_zscore_rmbatch_t = data.frame(t(training_zscore_rmbatch))
range(training_zscore_rmbatch_t) # -4.44812 4817.93667



platform_train = model_clinic[inTrain,'platform']


group_train = model_clinic[inTrain,'Canon']

tss_z_pCA <- p_cal(trainTransformed = training_zscore_rmbatch,
                   trainLabels = group_train$Canon,
                   trainLabels_sig = 'CA',paired = FALSE)# 2.30879442184768e-10   0.999946126323105  

dim(tss_z_pCA) # 43352     2
tss_z_pCAcut <- tss_z_pCA[which(tss_z_pCA[,1] < as.numeric(0.01) ),]
dim(tss_z_pCAcut)  #8527    2
### choose top 300
tss_z_pCAcut <- tss_z_pCAcut[order(tss_z_pCAcut$tout),]
tss_z_pCAcut <- tss_z_pCAcut[1:as.numeric(300),]

###
dim(tss_z_pCAcut) #  3239    2

### build the model based on features:tss_z_pCAcut$V2

train_topbins<-as.data.frame(training_zscore_rmbatch)[,tss_z_pCAcut$V2]
dim(train_topbins) #  339 3239
test_topbins<-as.data.frame(testing_zscore)[,tss_z_pCAcut$V2]
dim(test_topbins) # 112 3239

train_topbins$group <- factor(trainLabels$Canon,levels = c('CA','nonCA'))
test_topbins$group <- factor(testLabels$Canon,levels = c('CA','nonCA'))



### build model based on rf
set.seed(123) 
ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 3,
                    allowParallel = TRUE,classProbs = TRUE,savePredictions = TRUE,search='grid',summaryFunction = twoClassSummary)
set.seed(123) 
ctrl$sampling <- "up"
tunegrid <- expand.grid(mtry = c(1))
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

modeltss_CAnonCA <- rf_default1
#tss_CAnonCA1000 <- rf_default1
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


test_obs <- data.frame(test_topbins$group)
test_obs$obs_prob <- 'm'
test_obs[which(test_obs$test_topbins.group=='nonCA'),]$obs_prob = 0
test_obs[which(test_obs$test_topbins.group=='CA'),]$obs_prob = 1
roc_data <- roc(test_obs$obs_prob, pred[,1],type="prob")
print(roc_data)
#print(c(args[1],args[2],args[3],args[4]))
save(alltss_zscore,tss_z_dhpbatchcut,DH_zscore_rmbatch_t,
     training_zscore_rmbatch,tss_z_pCAcut,train_topbins,test_topbins,
     modeltss_CAnonCA,pred,roc_data,
     file='/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/tssCA.Rdata')
####################################################################################
############################         H vs  inflammation   ##########################
####################################################################################
load('/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/tssCA.Rdata')
### use the data rm batch (DH_zscore_rmbatch_t) clinic (model_clinic)
tssrmbatch <- DH_zscore_rmbatch_t
model_clinic <- read_excel("~/Desktop/XJ/XJ_paper/data/clinic/model_clinic2.xlsx")
###only WI and H included
model_clinic_WIH <- model_clinic[which(model_clinic$disease !='T'),]
tssrmbatch <- tssrmbatch[,model_clinic_WIH$sample]
tssrmbatch <- data.frame(t(tssrmbatch))
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
tss_z_pCAcut <- tss_z_pCA[which(tss_z_pCA[,1] < as.numeric(0.01) ),]
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
modeltss_WIH <- train(group~., 
                     data=train_topbins, 
                     method='rf', 
                     #metric='Accuracy',  #Metric compare model is Accuracy
                     metric = "ROC",
                     tuneGrid=tunegrid,
                     ntree=c(500),
                     trControl=ctrl)

rf_default1 <- modeltss_WIH
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


pred_tssWIH <- pred
roc_data_tssWIH <- roc_data
train_topbins_tssWIH <- train_topbins 
test_topbins_tssWIH<- test_topbins


save(train_topbins_tssWIH,test_topbins_tssWIH,
     modeltss_WIH,pred_tssWIH,roc_data_tssWIH,
     file='/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/tssWIH.Rdata')


####################################################################################
############################    localization: G vs  C ##############################
####################################################################################
load('/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/tssCA.Rdata')
### use the data rm batch (DH_zscore_rmbatch_t) clinic (model_clinic)
tssrmbatch <- DH_zscore_rmbatch_t
model_clinic <- read_excel("~/Desktop/XJ/XJ_paper/data/clinic/model_clinic2.xlsx")
###only WI and H included
model_clinic_WC <- model_clinic[which(model_clinic$disease =='T'),]
tssrmbatch <- tssrmbatch[,model_clinic_WC$sample]
tssrmbatch <- data.frame(t(tssrmbatch))
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
tss_z_pCAcut <- tss_z_pCA[which(tss_z_pCA[,1] < as.numeric(0.01) ),]
dim(tss_z_pCAcut)  # 1201    2
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
modeltss_WC <- train(group~., 
                      data=train_topbins, 
                      method='rf', 
                      #metric='Accuracy',  #Metric compare model is Accuracy
                      metric = "ROC",
                      tuneGrid=tunegrid,
                      ntree=c(500),
                      trControl=ctrl)

rf_default1 <- modeltss_WC
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


pred_tssWC <- pred
roc_data_tssWC <- roc_data
train_topbins_tssWC <- train_topbins 
test_topbins_tssWC<- test_topbins


save(train_topbins_tssWC,test_topbins_tssWC,
     modeltss_WC,pred_tssWC,roc_data_tssWC,
     file='/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/tssWC.Rdata')





#cnv_CAnonpred,cnv_roc_data,cnv_test_topbins,cnv_train_topbins,
#model_clinic3_cnv,model_cnvCAnonCA,model_data_cnv250gl_t,













