# Rscript xxx.R ''
args=commandArgs(T)
### import library
# Rscript /home/pingyi/0913/XJ_succ/model_bash_run_CAnonCA.R '~/0913/model_data_cnv250gl.csv' '~/0913/model_clinic.csv' 
# (0.01 0.5 250 +3)
# args[1] is output Rdata '~/0913/model_output_Rdata/cnv_model_CAnonCA.Rdata'
# args[2] is input data set samples are col 450 do not include 'pos'
# eg. '~/0913/model_data_cnv250gl.csv'
# args[3] is the clinic dataset eg. '~/0913/model_clinic.csv'
# args[4] is colname of group in clinic eg.'Canon'
# args[5] is positive of group eg. 'CA'

print(paste0(args[1],args[2],args[3]))

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
model_data_cnv250gl <- model_data_cnv250gl[,model_clinic3_cnv$sample]

### Training && test
model_data_cnv250gl_t <- data.frame(t(model_data_cnv250gl))
### 75% for training 25% for testing per_pa_rmbatch
set.seed(123)
inTrain=createDataPartition(y=model_clinic3_cnv$Canon,p=as.numeric(0.75),list=FALSE)
training=model_data_cnv250gl_t[inTrain,]
dim(training) # 338 3822
testing=model_data_cnv250gl_t[-inTrain,]
dim(testing)
### save the group of training && testing 
trainLabels <- data.frame(Canon=model_clinic3_cnv[inTrain,'Canon'])
testLabels <- data.frame(Canon=model_clinic3_cnv[-inTrain,'Canon'])
test_clinic_matrix<- model_clinic3_cnv[-inTrain,]
dim(training) #  338 90717
dim(testing)  # 112 90717

###

train_topbins<-as.data.frame(training)
dim(train_topbins) #  339 3239
test_topbins<-as.data.frame(testing)
dim(test_topbins) # 112 3239


train_topbins$group <- factor(trainLabels$Canon,levels = c('CA','nonCA'))
train_topbins$group <- factor(testLabels$Canon,levels = c('CA','nonCA'))


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
model_cnvCAnonCA <- rf_default1
print(rf_default1)

pred <- predict(rf_default1,newdata = test_topbins,type = "prob")
pred$sample = rownames(pred)
model_clinic = model_clinic3_cnv
pred <- left_join(pred,model_clinic,by='sample')
pred$res = 1
pred[which(pred$CA<0.5),'res'] = 0 

pred$true = 1
pred[which(pred$Canon != 'CA'),'true'] = 0 

pred$equ = 1
pred[which(pred$res != pred$true),'equ'] = 0 

test_obs <- data.frame(test_topbins$group)
test_obs$obs_prob <- 'm'
test_obs[which(test_obs$test_topbins.group=='nonCA'),]$obs_prob = 0
test_obs[which(test_obs$test_topbins.group=='CA'),]$obs_prob = 1
roc_data <- roc(test_obs$obs_prob, pred[,1],type="prob")

roc_data


######################################################################
############################### Save data ############################
######################################################################

save(model_data_cnv250gl_t,model_clinic3_cnv,
     train_topbins,test_topbins,model_cnvCAnonCA,pred,roc_data,
     file = args[1])


######################################################################






