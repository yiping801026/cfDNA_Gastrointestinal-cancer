# Rscript xxx.R ''
args=commandArgs(T)
### import library
### the final one :
# nohup Rscript /home/pingyi/0913/XJ_succ/model_bash_run_WIH.R '~/0913/cnvWIH.Rdata' 
# '~/0913/CNV_rmb_pca.csv' '~/0913/model_clinic.csv' 0.05 250 40 > report40.log 2>&1
# args[1] is output Rdata '~/0913/model_output_Rdata/cnv_model_CAnonCA.Rdata'
# args[2] is input data set samples are col 450 do not include 'pos'
# eg. '~/0913/model_data_cnv250gl.csv'
# args[3] is the clinic dataset eg. '~/0913/model_clinic.csv'
# args[4] is 0.05 for group p cutoff
# args[5] is the number of bins included in the model 250
# args[6] is the three gainloss median var included in the model eg."~/0913/count_100.csv"

### import inhouse functions
source('/home/pingyi/XJ_model/bash_run/script/used_function.R')
#args[1] = '/home/pingyi/0913/XJ_succ/model_bash_run_WIH.R'

print(paste0(args[1],args[2],args[3],args[4],args[5],args[6]))

library(readxl)
library(dplyr)
library(factoextra) 
library(ggplot2)
library(pROC)
library(stats)
library(caret)
######################################################################
############################### Run model ############################
######################################################################
model_data_cnv250gl =  read.csv(args[2],row.names = 1)
model_clinic3_cnv = read.csv(args[3],row.names = 1)
cnv_gainloss <- read.csv("~/0913/count_100.csv",row.names = 1)
model_clinic3_cnv <- model_clinic3_cnv[model_clinic3_cnv$disease  !='T', ]
model_data_cnv250gl <- model_data_cnv250gl[,model_clinic3_cnv$sample]

### Training && test
model_data_cnv250gl_t <- data.frame(t(model_data_cnv250gl))
### 75% for training 25% for testing per_pa_rmbatch
set.seed(123)
inTrain=createDataPartition(y=model_clinic3_cnv$disease,p=as.numeric(0.75),list=FALSE)
training=model_data_cnv250gl_t[inTrain,]
### use median & variance & percentage


dim(training) # 338 3822
testing=model_data_cnv250gl_t[-inTrain,]
dim(testing)
### save the group of training && testing 
trainLabels <- data.frame(disease=model_clinic3_cnv[inTrain,'disease'])
testLabels <- data.frame(disease=model_clinic3_cnv[-inTrain,'disease'])
test_clinic_matrix<- model_clinic3_cnv[-inTrain,]
dim(training) #  338 90717
dim(testing)  # 112 90717
##### t-test for CAnon

group_train = data.frame(disease = model_clinic3_cnv[inTrain,'disease'])
training_t <- data.frame(t(training))
cnv_z_pCA <- p_cal(trainTransformed = training,
                   trainLabels = group_train$disease,
                   trainLabels_sig = 'I',paired = FALSE)# 2.30879442184768e-10   0.999946126323105  

cnv_z_pCAcut <- cnv_z_pCA[which(cnv_z_pCA[,1] < as.numeric(args[4]) ),]
dim(cnv_z_pCAcut) # 3457    2

### choose top bins
tss_z_pCAcut <- cnv_z_pCAcut
tss_z_pCAcut <- tss_z_pCAcut[order(tss_z_pCAcut$tout),]
tss_z_pCAcut <- tss_z_pCAcut[1: as.numeric(args[5]),]


train_topbins = data.frame(training)[,tss_z_pCAcut$V2]
add_fea_train = cnv_gainloss[inTrain,]
train_topbins <- cbind(train_topbins,add_fea_train[,c(1,3,12)])
dim(train_topbins) #  197 253

test_topbins = data.frame(testing)[,tss_z_pCAcut$V2]
add_fea_test = cnv_gainloss[which(cnv_gainloss$sample %in% rownames(test_topbins)),]
test_topbins <- cbind(test_topbins,add_fea_test[,c(1,3,12)])
dim(test_topbins) # 65

train_topbins$group <- factor(trainLabels$disease,levels = c('I','H'))
test_topbins$group <- factor(testLabels$disease,levels = c('I','H'))

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
                     ntree=c(1000),
                     trControl=ctrl)
rf_default1$bestTune
model_cnvWIH <- rf_default1
print(rf_default1)

pred <- predict(rf_default1,newdata = test_topbins,type = "prob")
pred$sample = rownames(pred)
model_clinic = model_clinic3_cnv
pred <- left_join(pred,model_clinic,by='sample')
pred$res = 1
pred[which(pred$I<0.5),'res'] = 0 

pred$true = 1
pred[which(pred$disease != 'I'),'true'] = 0 

pred$equ = 1
pred[which(pred$res != pred$true),'equ'] = 0 

test_obs <- data.frame(test_topbins$group)
test_obs$obs_prob <- 'm'
test_obs[which(test_obs$test_topbins.group=='H'),]$obs_prob = 0
test_obs[which(test_obs$test_topbins.group=='I'),]$obs_prob = 1
roc_data <- roc(test_obs$obs_prob, pred[,1],type="prob")

roc_data


######################################################################
############################### Save data ############################
######################################################################

save(model_data_cnv250gl_t,model_clinic3_cnv,
     train_topbins,test_topbins,model_cnvWIH,pred,roc_data,
     file = args[1])


######################################################################






