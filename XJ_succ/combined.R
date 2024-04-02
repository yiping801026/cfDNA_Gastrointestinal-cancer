library(dplyr)
library(pROC)
library(ggplot2)
library(dplyr)

load(file = "/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/cnv_model_CAnonCA.Rdata")
load(file = "/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/cnv_model_WC.Rdata")
load(file = "/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/cnv_model_WIH41.Rdata")
load(file = '/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/fragCA.Rdata')
load(file = '/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/fragWIH.Rdata')
load(file = '/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/fragWC.Rdata')
load(file = '/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/tssCA.Rdata')
load(file = '/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/tssWC.Rdata')
load(file = '/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/tssWIH.Rdata')

WC_test
frag_CAnopred
cnv_CAnonpred 
cnv_CAnonpred <- pred##tss
cnv_roc_data <- roc_data
cnv_test_topbins<- test_topbins
cnv_train_topbins<- train_topbins

save(cnv_CAnonpred,cnv_roc_data,cnv_test_topbins,cnv_train_topbins,
     model_clinic3_cnv,model_cnvCAnonCA,model_data_cnv250gl_t,
     file = "/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/cnv_model_CAnonCA.Rdata")


cnv_WC_pred <- pred
cnv_WCroc_data <- roc_data
cnvWC_test_topbins<- test_topbins
cnvWC_train_topbins<- train_topbins
save(model_clinic3_cnv,cnv_WC_pred,model_cnvWC,model_data_cnv250gl_t,cnv_WCroc_data,
  cnvWC_test_topbins,cnvWC_train_topbins,file="/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/cnv_model_WC.Rdata")


cnv_WIH_pred <- pred
cnv_WIHroc_data <- roc_data
cnvWIH_test_topbins<- test_topbins
cnvWIH_train_topbins<- train_topbins
save(model_clinic3_cnv,cnv_WIH_pred,model_cnvWIH,model_data_cnv250gl_t,cnv_WIHroc_data,
     cnvWIH_test_topbins,cnvWIH_train_topbins,file = "/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/cnv_model_WIH41.Rdata")


####################################################################################
############################        CA nonCA              ##########################
####################################################################################
###### trining cohort 
### training results of several models
tss_CAnon_trare <- data.frame(modeltss_CAnonCA$pred)
tss_CAnon_trare <- tss_CAnon_trare[which(tss_CAnon_trare$mtry == 1),]
frag_CAnon_trare <- data.frame(modelfrag_CAnonCA$pred) ### mtry : 20-25, so choose the mtry == 25
frag_CAnon_trare <- frag_CAnon_trare[which(frag_CAnon_trare$mtry == 5),]
cnv_CAnon_trare <- data.frame(model_cnvCAnonCA $pred)  ### mtry : 20-25, so choose the mtry == 25
cnv_CAnon_trare <- cnv_CAnon_trare[which(cnv_CAnon_trare$mtry == 3),]



modeltss_CAnonCA # 0.9175821  0.7042857  0.8988596
modelfrag_CAnonCA # 5     0.7953425  0.6122222  0.8336842
model_cnvCAnonCA # 3     0.8600524  0.6998413  0.8493860
### combined model
#pred_tssWC

### for trianing data -- mean of 3 repeats
roc_mean_calculate <- function(matrix){
  ### ROC result of the original models
  matrix$obs_num <- rep(1,nrow(matrix))
  matrix[matrix$obs == 'nonCA','obs_num'] = 0
  matrix_result <- aggregate(matrix$CA, by = list(matrix$rowIndex), FUN = mean)
  matrix_resgroup <- unique(data.frame(rowIndex=matrix$rowIndex,obs_num=matrix$obs_num))
  colnames(matrix_result)= c('rowIndex','CA')
  matrix_result <- left_join(matrix_result,matrix_resgroup,by='rowIndex')
  roc_data <- roc(matrix_result$obs_num, matrix_result$CA,type="prob")
  return(roc_data)
}
### for training data -- mean of 3 repeats
roc_mean_calculate_matrix <- function(matrix){
  ### ROC result of the original models
  matrix$obs_num <- rep(1,nrow(matrix))
  matrix[matrix$obs == 'nonCA','obs_num'] = 0
  matrix_result <- aggregate(matrix$CA, by = list(matrix$rowIndex), FUN = mean)
  matrix_resgroup <- unique(data.frame(rowIndex=matrix$rowIndex,obs_num=matrix$obs_num))
  colnames(matrix_result)= c('rowIndex','CA')
  matrix_result <- left_join(matrix_result,matrix_resgroup,by='rowIndex')
  roc_data <- roc(matrix_result$obs_num, matrix_result$CA,type="prob")
  return(matrix_result)
}
### roc for models
cnv_roc_result <- roc_mean_calculate(cnv_CAnon_trare)
tss_roc_result <- roc_mean_calculate(tss_CAnon_trare)
frag_roc_result <- roc_mean_calculate(frag_CAnon_trare)

### pred matrix for models
cnv_rocplot <- roc_mean_calculate_matrix(cnv_CAnon_trare)
tss_rocplot <- roc_mean_calculate_matrix(tss_CAnon_trare)
frag_rocplot <- roc_mean_calculate_matrix(frag_CAnon_trare)

### for combined model
combined_roc <- left_join(cnv_rocplot,tss_rocplot,by='rowIndex')
combined_roc <- left_join(combined_roc,frag_rocplot,by='rowIndex')
combined_roc <- data.frame(combined_roc)
combined_roc$mean <- apply(data.frame(combined_roc[,c(2,4,6)]), 1, mean)
combinde_auc <- roc(combined_roc$obs_num, combined_roc$mean,type="prob")
rownames(combined_roc) <- rownames(training)
write.csv(combined_roc,'/Users/pingyi/Desktop/XJ_new/script/out_matrix/fig5.1CAnonCAtraining_roc.csv')

### four models together
combined_roc <- combined_roc[,c(1,2,4,6,7,8)]
colnames(combined_roc)[c(2:4,6)] <- c('cnv','tss','frag','combined')
### make the multi-plot of ROC curves for CA vs nonCA of 4 models -- training dataset
# 0.8600524 0.9175821 0.7953425 0.886
roc.list_train <- roc(obs_num ~ cnv + tss + frag + combined, data = combined_roc,type="prob")
g.list <- ggroc(roc.list_train,size = 1,legacy.axes = TRUE) ### change the line size
p <- g.list + scale_colour_manual(values = c("#2878b5", "#9ac9db", "#f8ac8c",'#c82423'), ### change line colors
                             labels = c('cnv (0.860)','tss (0.918)','frag (0.795)','combined (0.886)'))+ ### change legend names
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")+ ### add middle line
  labs(title = 'Training Cohort',col = 'Models (AUC)', x='False Positive Rate',y='True Positive Rate')+ ### change legend title(col) & x,ynames 
  ggplot2::theme(
    panel.grid.major = element_line(colour = NA),
    panel.background = element_rect(fill = "transparent", colour = 'black'),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.minor = element_blank()
  ) 

p + ggtitle("Screening Model (CA vs nonCA)") +
  ggplot2::theme_bw() +
  #ggplot2::ylim(-ylimn,ylimn)+
  #ggplot2::xlim(-xlimn,xlimn)+
  coord_fixed(ratio = 1)


############################# test


###### test cohort 
### test results of several model
frag_CAnopred
pred##tss
cnv_CAnonpred
#CAnonCA_testroc$try <- cnv_CAnonpred$disease
CAnonCA_testroc<- data.frame(sample_names=rownames(cnv_CAnonpred),obs_num = c(rep(1,47),rep(0,65)),
                             cnv = cnv_CAnonpred[,1],tss=pred[,1],frag=frag_CAnopred[,1])

CAnonCA_testroc$combined <- apply(data.frame(CAnonCA_testroc[,3:5]), 1, mean)
combinde_auc <- roc(CAnonCA_testroc$obs_num, CAnonCA_testroc$combined,type="prob")

write.csv(CAnonCA_testroc,'/Users/pingyi/Desktop/XJ_new/script/out_matrix/fig5.2CAnonCA_testroc.csv')
### calculate the auc for 4 models

#test_cnv <- roc(CAnonCA_testroc$obs_num,CAnonCA_testroc$cnv,type="prob") #0.9579
#test_tss <- roc(CAnonCA_testroc$obs_num,CAnonCA_testroc$tss,type="prob") #1
#test_frag <- roc(CAnonCA_testroc$obs_num,CAnonCA_testroc$frag,type="prob") #0.8582
#test_combined <- roc(CAnonCA_testroc$obs_num,CAnonCA_testroc$combined,type="prob") #1
 0.8653
### make the multi-plot of ROC curves for CA vs nonCA of 4 models -- training dataset
roc.list_test <- roc(obs_num ~ cnv + tss + frag + combined, data = CAnonCA_testroc,type="prob")
g.list2 <- ggroc(roc.list_test,size = 1,legacy.axes = TRUE) ### change the line size
g.list2 + scale_colour_manual(values = c("#2878b5", "#9ac9db", "#f8ac8c",'#c82423'), ### change line colors
                              labels = c('cnv (0.925)','tss (0.937)','frag (0.865)','combined (0.9381)'))+ ### change legend names
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")+ ### add middle line
  labs(title = 'Validation Cohort',col = 'Models (AUC)', x='False Positive Rate',y='True Positive Rate')+ ### change legend title(col) & x,ynames 
  ggplot2::theme(
    panel.grid.major = element_line(colour = NA),
    panel.background = element_rect(fill = "transparent", colour = 'black'),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.minor = element_blank()
  )+ ggtitle("Screening Model (CA vs nonCA)") +
  ggplot2::theme_bw() +
  #ggplot2::ylim(-ylimn,ylimn)+
  #ggplot2::xlim(-xlimn,xlimn)+
  coord_fixed(ratio = 1)

####################################################################################
############################        H vs  inflammation    ##########################
####################################################################################
#tss 27

tss_CAnon_trare <- data.frame(modeltss_WIH$pred)
tss_CAnon_trare <- tss_CAnon_trare[which(tss_CAnon_trare$mtry == 3),]
frag_CAnon_trare <- data.frame(modelfrag_WIH$pred) ### mtry : 20-25, so choose the mtry == 25
frag_CAnon_trare <- frag_CAnon_trare[which(frag_CAnon_trare$mtry == 2),]
cnv_CAnon_trare <- data.frame(model_cnvWIH $pred)  ### mtry : 20-25, so choose the mtry == 25
cnv_CAnon_trare <- cnv_CAnon_trare[which(cnv_CAnon_trare$mtry == 4),]

### trining cohort 
### training results of several models
model_cnvWIH
modelfrag_WIH
modeltss_WIH
cnv 0.8801435 tss 0.8856074 frag 0.8641636 c 0.9513

### for trining data -- mean of 3 repeats
roc_mean_calculate <- function(matrix){
  ### ROC result of the original models
  #matrix = tss_CAnon_trare
  matrix$obs_num <- rep(1,nrow(matrix))
  matrix[matrix$obs == 'H','obs_num'] = 0
  matrix_result <- aggregate(matrix$I, by = list(matrix$rowIndex), FUN = mean)
  matrix_resgroup <- unique(data.frame(rowIndex=matrix$rowIndex,obs_num=matrix$obs_num))
  colnames(matrix_result)= c('rowIndex','I')
  matrix_result <- left_join(matrix_result,matrix_resgroup,by='rowIndex')
  roc_data <- roc(matrix_result$obs_num, matrix_result$I,type="prob")
  return(roc_data)
}
### for trining data -- mean of 3 repeats
roc_mean_calculate_matrix <- function(matrix){
  #matrix=tss_CAnon_trare
  ### ROC result of the original models
  matrix$obs_num <- rep(1,nrow(matrix))
  matrix[matrix$obs == 'H','obs_num'] = 0
  matrix_result <- aggregate(matrix$I, by = list(matrix$rowIndex), FUN = mean)
  matrix_resgroup <- unique(data.frame(rowIndex=matrix$rowIndex,obs_num=matrix$obs_num))
  colnames(matrix_result)= c('rowIndex','I')
  matrix_result <- left_join(matrix_result,matrix_resgroup,by='rowIndex')
  roc_data <- roc(matrix_result$obs_num, matrix_result$I,type="prob")
  return(matrix_result)
}
### roc for models
tss_roc_result <- roc_mean_calculate(tss_CAnon_trare)
frag_roc_result <- roc_mean_calculate(frag_CAnon_trare)
cnv_roc_result <- roc_mean_calculate(cnv_CAnon_trare)

### pred matrix for models
tss_rocplot <- roc_mean_calculate_matrix(tss_CAnon_trare)
frag_rocplot <- roc_mean_calculate_matrix(frag_CAnon_trare)
cnv_rocplot <- roc_mean_calculate_matrix(cnv_CAnon_trare)

### for combined model
combined_roc <- left_join(tss_rocplot,frag_rocplot,by='rowIndex')
combined_roc <- left_join(combined_roc,cnv_rocplot,by='rowIndex')

combined_roc <- data.frame(combined_roc[,c(1:2,4,6:7)])
colnames(combined_roc) <- c("rowIndex" ,'tss','frag','cnv','obs_num')
combined_roc$mean <- apply(data.frame(combined_roc[,c(2,3,4)]), 1, mean)
combinde_auc <- roc(combined_roc$obs_num, combined_roc$mean,type="prob")
rownames(combined_roc) <- rownames(train_topbins_fragWIH)
write.csv(combined_roc,'/Users/pingyi/Desktop/XJ_new/script/out_matrix/fig5.1WIH_roc.csv')
colnames(combined_roc) <- c("rowIndex" ,'tss','frag','cnv','obs_num','combined')
### make the multi-plot of ROC curves for CA vs nonCA of 4 models -- training dataset

roc.list_train <- roc(obs_num ~ tss + frag+cnv + combined, data = combined_roc,type="prob")
g.list <- ggroc(roc.list_train,size = 1,legacy.axes = TRUE) ### change the line size
g.list + scale_colour_manual(values = c( "#2878b5", "#9ac9db", "#f8ac8c",'#c82423'), ### change line colors
                             labels = c('cnv (0.880)','tss (0.886)','frag (0.864)','combined (0.951)'))+ ### change legend names
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")+ ### add middle line
  labs(title = 'Training Cohort',col = 'Models (AUC)', x='False Positive Rate',y='True Positive Rate')+ ### change legend title(col) & x,ynames 
  ggplot2::theme(
    panel.grid.major = element_line(colour = NA),
    panel.background = element_rect(fill = "transparent", colour = 'black'),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.minor = element_blank()
  ) + ggtitle("Screening Model (Inflam vs Healthy)") +
  ggplot2::theme_bw() +
  #ggplot2::ylim(-ylimn,ylimn)+
  #ggplot2::xlim(-xlimn,xlimn)+
  coord_fixed(ratio = 1)

################ test cohort
cnv_WIH_pred
pred_fragWIH
pred_tssWIH
### test results of several models
WIH_testroc<- data.frame(sample_names=rownames(cnv_WIH_pred),obs_num = c(rep(1,23),rep(0,42)),
                         tss=pred_tssWIH[,1],frag=pred_fragWIH[,1],cnv= cnv_WIH_pred[,1])

WIH_testroc$combined <- apply(data.frame(WIH_testroc[,3:5]), 1, mean)
write.csv(WIH_testroc,'/Users/pingyi/Desktop/XJ_new/script/out_matrix/fig5.2WIH_testroc.csv')
#WIH_testroc$try = pred_tssWIH$disease
### calculate the auc for 4 models

### make the multi-plot of ROC curves for CA vs nonCA of 4 models -- training dataset
roc.list_test <- roc(obs_num ~ tss + frag+cnv + combined, data = WIH_testroc,type="prob")
g.list2 <- ggroc(roc.list_test,size = 1,legacy.axes = TRUE) ### change the line size
g.list2 + scale_colour_manual(values = c("#2878b5", "#9ac9db", "#f8ac8c",'#c82423'), ### change line colors
                              labels = c('cnv (0.851)','tss (0.798)','frag (0.755)','combined (0.890)'))+ ### change legend names
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")+ ### add middle line
  labs(title = 'Validation Cohort',col = 'Models (AUC)', x='False Positive Rate',y='True Positive Rate')+ ### change legend title(col) & x,ynames 
  ggplot2::theme(
    panel.grid.major = element_line(colour = NA),
    panel.background = element_rect(fill = "transparent", colour = 'black'),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.minor = element_blank()
  ) + ggtitle("Screening Model (Inflam vs Healthy)") +
  ggplot2::theme_bw() +
  #ggplot2::ylim(-ylimn,ylimn)+
  #ggplot2::xlim(-xlimn,xlimn)+
  coord_fixed(ratio = 1)


####################################################################################
############################        CA nonCA              ##########################
####################################################################################
###### trining cohort 
### training results of several models
tss_CAnon_trare <- data.frame(modeltss_CAnonCA$pred)
tss_CAnon_trare <- tss_CAnon_trare[which(tss_CAnon_trare$mtry == 1),]
frag_CAnon_trare <- data.frame(modelfrag_CAnonCA$pred) ### mtry : 20-25, so choose the mtry == 25
frag_CAnon_trare <- frag_CAnon_trare[which(frag_CAnon_trare$mtry == 5),]
cnv_CAnon_trare <- data.frame(model_cnvCAnonCA $pred)  ### mtry : 20-25, so choose the mtry == 25
cnv_CAnon_trare <- cnv_CAnon_trare[which(cnv_CAnon_trare$mtry == 3),]



modeltss_CAnonCA # 0.9175821  0.7042857  0.8988596
modelfrag_CAnonCA # 5     0.7953425  0.6122222  0.8336842
model_cnvCAnonCA # 3     0.8600524  0.6998413  0.8493860
### combined model
#pred_tssWC

### for trianing data -- mean of 3 repeats
roc_mean_calculate <- function(matrix){
  ### ROC result of the original models
  matrix$obs_num <- rep(1,nrow(matrix))
  matrix[matrix$obs == 'nonCA','obs_num'] = 0
  matrix_result <- aggregate(matrix$CA, by = list(matrix$rowIndex), FUN = mean)
  matrix_resgroup <- unique(data.frame(rowIndex=matrix$rowIndex,obs_num=matrix$obs_num))
  colnames(matrix_result)= c('rowIndex','CA')
  matrix_result <- left_join(matrix_result,matrix_resgroup,by='rowIndex')
  roc_data <- roc(matrix_result$obs_num, matrix_result$CA,type="prob")
  return(roc_data)
}
### for training data -- mean of 3 repeats
roc_mean_calculate_matrix <- function(matrix){
  ### ROC result of the original models
  matrix$obs_num <- rep(1,nrow(matrix))
  matrix[matrix$obs == 'nonCA','obs_num'] = 0
  matrix_result <- aggregate(matrix$CA, by = list(matrix$rowIndex), FUN = mean)
  matrix_resgroup <- unique(data.frame(rowIndex=matrix$rowIndex,obs_num=matrix$obs_num))
  colnames(matrix_result)= c('rowIndex','CA')
  matrix_result <- left_join(matrix_result,matrix_resgroup,by='rowIndex')
  roc_data <- roc(matrix_result$obs_num, matrix_result$CA,type="prob")
  return(matrix_result)
}
### roc for models
cnv_roc_result <- roc_mean_calculate(cnv_CAnon_trare)
tss_roc_result <- roc_mean_calculate(tss_CAnon_trare)
frag_roc_result <- roc_mean_calculate(frag_CAnon_trare)

### pred matrix for models
cnv_rocplot <- roc_mean_calculate_matrix(cnv_CAnon_trare)
tss_rocplot <- roc_mean_calculate_matrix(tss_CAnon_trare)
frag_rocplot <- roc_mean_calculate_matrix(frag_CAnon_trare)

### for combined model
combined_roc <- left_join(cnv_rocplot,tss_rocplot,by='rowIndex')
combined_roc <- left_join(combined_roc,frag_rocplot,by='rowIndex')
combined_roc <- data.frame(combined_roc)
combined_roc$mean <- apply(data.frame(combined_roc[,c(2,4,6)]), 1, mean)
combinde_auc <- roc(combined_roc$obs_num, combined_roc$mean,type="prob")
rownames(combined_roc) <- rownames(training)
write.csv(combined_roc,'/Users/pingyi/Desktop/XJ_new/script/out_matrix/fig5.1CAnonCAtraining_roc.csv')

### four models together
combined_roc <- combined_roc[,c(1,2,4,6,7,8)]
colnames(combined_roc)[c(2:4,6)] <- c('cnv','tss','frag','combined')
### make the multi-plot of ROC curves for CA vs nonCA of 4 models -- training dataset
# 0.8600524 0.9175821 0.7953425 0.886
roc.list_train <- roc(obs_num ~ cnv + tss + frag + combined, data = combined_roc,type="prob")
g.list <- ggroc(roc.list_train,size = 1,legacy.axes = TRUE) ### change the line size
p <- g.list + scale_colour_manual(values = c("#2878b5", "#9ac9db", "#f8ac8c",'#c82423'), ### change line colors
                                  labels = c('cnv (0.860)','tss (0.918)','frag (0.795)','combined (0.886)'))+ ### change legend names
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")+ ### add middle line
  labs(title = 'Training Cohort',col = 'Models (AUC)', x='False Positive Rate',y='True Positive Rate')+ ### change legend title(col) & x,ynames 
  ggplot2::theme(
    panel.grid.major = element_line(colour = NA),
    panel.background = element_rect(fill = "transparent", colour = 'black'),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.minor = element_blank()
  ) 

p + ggtitle("Screening Model (CA vs nonCA)") +
  ggplot2::theme_bw() +
  #ggplot2::ylim(-ylimn,ylimn)+
  #ggplot2::xlim(-xlimn,xlimn)+
  coord_fixed(ratio = 1)


############################# test


###### test cohort 
### test results of several model
frag_CAnopred


pred##tss
cnv_CAnonpred
#CAnonCA_testroc$try <- cnv_CAnonpred$disease
CAnonCA_testroc<- data.frame(sample_names=rownames(cnv_CAnonpred),obs_num = c(rep(1,47),rep(0,65)),
                             cnv = cnv_CAnonpred[,1],tss=pred[,1],frag=frag_CAnopred[,1])

CAnonCA_testroc$combined <- apply(data.frame(CAnonCA_testroc[,3:5]), 1, mean)
combinde_auc <- roc(CAnonCA_testroc$obs_num, CAnonCA_testroc$combined,type="prob")

write.csv(CAnonCA_testroc,'/Users/pingyi/Desktop/XJ_new/script/out_matrix/fig5.2CAnonCA_testroc.csv')
### calculate the auc for 4 models

#test_cnv <- roc(CAnonCA_testroc$obs_num,CAnonCA_testroc$cnv,type="prob") #0.9579
#test_tss <- roc(CAnonCA_testroc$obs_num,CAnonCA_testroc$tss,type="prob") #1
#test_frag <- roc(CAnonCA_testroc$obs_num,CAnonCA_testroc$frag,type="prob") #0.8582
#test_combined <- roc(CAnonCA_testroc$obs_num,CAnonCA_testroc$combined,type="prob") #1
0.8653
### make the multi-plot of ROC curves for CA vs nonCA of 4 models -- training dataset
roc.list_test <- roc(obs_num ~ cnv + tss + frag + combined, data = CAnonCA_testroc,type="prob")
g.list2 <- ggroc(roc.list_test,size = 1,legacy.axes = TRUE) ### change the line size
g.list2 + scale_colour_manual(values = c("#2878b5", "#9ac9db", "#f8ac8c",'#c82423'), ### change line colors
                              labels = c('cnv (0.925)','tss (0.937)','frag (0.865)','combined (0.9381)'))+ ### change legend names
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")+ ### add middle line
  labs(title = 'Validation Cohort',col = 'Models (AUC)', x='False Positive Rate',y='True Positive Rate')+ ### change legend title(col) & x,ynames 
  ggplot2::theme(
    panel.grid.major = element_line(colour = NA),
    panel.background = element_rect(fill = "transparent", colour = 'black'),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.minor = element_blank()
  )+ ggtitle("Screening Model (CA vs nonCA)") +
  ggplot2::theme_bw() +
  #ggplot2::ylim(-ylimn,ylimn)+
  #ggplot2::xlim(-xlimn,xlimn)+
  coord_fixed(ratio = 1)

####################################################################################
############################        W vs C    ##########################
####################################################################################
#tss 27

tss_CAnon_trare <- data.frame(modeltss_WC$pred)
tss_CAnon_trare <- tss_CAnon_trare[which(tss_CAnon_trare$mtry == 1),]
frag_CAnon_trare <- data.frame(modelfrag_WC$pred) ### mtry : 20-25, so choose the mtry == 25
frag_CAnon_trare <- frag_CAnon_trare[which(frag_CAnon_trare$mtry == 1),]
cnv_CAnon_trare <- data.frame(model_cnvWC$pred)  ### mtry : 20-25, so choose the mtry == 25
cnv_CAnon_trare <- cnv_CAnon_trare[which(cnv_CAnon_trare$mtry == 4),]

### trining cohort 
### training results of several models
cnv 0.8801435 tss 0.8856074 frag 0.8641636 c 0.9513

### for trining data -- mean of 3 repeats
roc_mean_calculate <- function(matrix){
  ### ROC result of the original models
  #matrix = tss_CAnon_trare
  matrix$obs_num <- rep(1,nrow(matrix))
  matrix[matrix$obs == 'C','obs_num'] = 0
  matrix_result <- aggregate(matrix$W, by = list(matrix$rowIndex), FUN = mean)
  matrix_resgroup <- unique(data.frame(rowIndex=matrix$rowIndex,obs_num=matrix$obs_num))
  colnames(matrix_result)= c('rowIndex','W')
  matrix_result <- left_join(matrix_result,matrix_resgroup,by='rowIndex')
  roc_data <- roc(matrix_result$obs_num, matrix_result$W,type="prob")
  return(roc_data)
}
### for trining data -- mean of 3 repeats
roc_mean_calculate_matrix <- function(matrix){
  #matrix=tss_CAnon_trare
  ### ROC result of the original models
  matrix$obs_num <- rep(1,nrow(matrix))
  matrix[matrix$obs == 'C','obs_num'] = 0
  matrix_result <- aggregate(matrix$W, by = list(matrix$rowIndex), FUN = mean)
  matrix_resgroup <- unique(data.frame(rowIndex=matrix$rowIndex,obs_num=matrix$obs_num))
  colnames(matrix_result)= c('rowIndex','W')
  matrix_result <- left_join(matrix_result,matrix_resgroup,by='rowIndex')
  roc_data <- roc(matrix_result$obs_num, matrix_result$W,type="prob")
  return(matrix_result)
}

### roc for models
cnv_roc_result <- roc_mean_calculate(cnv_CAnon_trare)
tss_roc_result <- roc_mean_calculate(tss_CAnon_trare)
frag_roc_result <- roc_mean_calculate(frag_CAnon_trare)

### pred matrix for models
cnv_rocplot <- roc_mean_calculate_matrix(cnv_CAnon_trare)
tss_rocplot <- roc_mean_calculate_matrix(tss_CAnon_trare)
frag_rocplot <- roc_mean_calculate_matrix(frag_CAnon_trare)

### for combined model
combined_roc <- left_join(cnv_rocplot,tss_rocplot,by='rowIndex')
combined_roc <- left_join(combined_roc,frag_rocplot,by='rowIndex')
combined_roc <- data.frame(combined_roc)
combined_roc$mean <- apply(data.frame(combined_roc[,c(2,4,6)]), 1, mean)
combinde_auc <- roc(combined_roc$obs_num, combined_roc$mean,type="prob")
rownames(combined_roc) <- rownames(train_topbins_fragWC)
write.csv(combined_roc,'/Users/pingyi/Desktop/XJ_new/script/out_matrix/fig5.1WCtraining_roc.csv')

### four models together
combined_roc <- combined_roc[,c(1,2,4,6,7,8)]
colnames(combined_roc)[c(2:4,6)] <- c('cnv','tss','frag','combined')
### make the multi-plot of ROC curves for CA vs nonCA of 4 models -- training dataset
# 0.8600524 0.9175821 0.7953425 0.886
roc.list_train <- roc(obs_num ~ cnv + tss + frag + combined, data = combined_roc,type="prob")
g.list <- ggroc(roc.list_train,size = 1,legacy.axes = TRUE) ### change the line size
p <- g.list + scale_colour_manual(values = c("#2878b5", "#9ac9db", "#f8ac8c",'#c82423'), ### change line colors
                                  labels = c('cnv (0.705)','tss (0.967)','frag (0.720)','combined (0.850)'))+ ### change legend names
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")+ ### add middle line
  labs(title = 'Training Cohort',col = 'Models (AUC)', x='False Positive Rate',y='True Positive Rate')+ ### change legend title(col) & x,ynames 
  ggplot2::theme(
    panel.grid.major = element_line(colour = NA),
    panel.background = element_rect(fill = "transparent", colour = 'black'),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.minor = element_blank()
  ) 

p + ggtitle("Localization Model (Gastric vs Colorectal)") +
  ggplot2::theme_bw() +
  #ggplot2::ylim(-ylimn,ylimn)+
  #ggplot2::xlim(-xlimn,xlimn)+
  coord_fixed(ratio = 1)

############################# test


###### test cohort 
### test results of several model
pred_fragWC
pred_tssWC
cnv_WC_pred

#CAnonCA_testroc$try <- cnv_CAnonpred$disease
CAnonCA_testroc<- data.frame(sample_names=rownames(pred_fragWC),obs_num = c(rep(0,12),rep(1,35)),
                             cnv = cnv_WC_pred[,1],tss=pred_tssWC[,1],frag=pred_fragWC[,1])

CAnonCA_testroc$combined <- apply(data.frame(CAnonCA_testroc[,3:5]), 1, mean)
combinde_auc <- roc(CAnonCA_testroc$obs_num, CAnonCA_testroc$combined,type="prob")

write.csv(CAnonCA_testroc,'/Users/pingyi/Desktop/XJ_new/script/out_matrix/fig5.2WC_testroc.csv')
### calculate the auc for 4 models

#test_cnv <- roc(CAnonCA_testroc$obs_num,CAnonCA_testroc$cnv,type="prob") #0.9579
#test_tss <- roc(CAnonCA_testroc$obs_num,CAnonCA_testroc$tss,type="prob") #1
#test_frag <- roc(CAnonCA_testroc$obs_num,CAnonCA_testroc$frag,type="prob") #0.8582
#test_combined <- roc(CAnonCA_testroc$obs_num,CAnonCA_testroc$combined,type="prob") #1
0.8653
combined<- roc(CAnonCA_testroc$obs_num,CAnonCA_testroc$combined,type="prob")
### make the multi-plot of ROC curves for CA vs nonCA of 4 models -- training dataset
#roc.list_test <- roc(obs_num ~ cnv + tss + frag + combined, data = CAnonCA_testroc,type="prob")



roc.list_test <-list(cnv_WCroc_data,roc_data_tssWC,roc_data_fragWC,combined)
g.list2 <- ggroc(roc.list_test,size = 1,legacy.axes = TRUE) ### change the line size
g.list2 + scale_colour_manual(values = c("#2878b5", "#9ac9db", "#f8ac8c",'#c82423'), ### change line colors
                              labels = c('cnv (0.740)','tss (0.606)','frag (0.639)','combined (0.693)'))+ ### change legend names
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")+ ### add middle line
  labs(title = 'Validation Cohort',col = 'Models (AUC)', x='False Positive Rate',y='True Positive Rate')+ ### change legend title(col) & x,ynames 
  ggplot2::theme(
    panel.grid.major = element_line(colour = NA),
    panel.background = element_rect(fill = "transparent", colour = 'black'),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.minor = element_blank()
  )+ ggtitle("Localization Model (Gastric vs Colorectal)") +
  ggplot2::theme_bw() +
  #ggplot2::ylim(-ylimn,ylimn)+
  #ggplot2::xlim(-xlimn,xlimn)+
  coord_fixed(ratio = 1)





####################################################################################
############################        stage   && GC         ##########################
####################################################################################
################ training
########
combined_roc <- read.csv('/Users/pingyi/Desktop/XJ_new/script/out_matrix/fig5.1CAnonCAtraining_roc.csv')
CAnonCA_testroc <- read.csv('/Users/pingyi/Desktop/XJ_new/script/out_matrix/fig5.2CAnonCA_testroc.csv')
model_clinic <- read.csv("~/Desktop/XJ_succ/input_data/clinic/model_clinic.csv")

XJ_clinic_raw <- model_clinic
#XJ_clinic_raw <- read.csv('/Users/pingyi/Desktop/XJ_new/script/out_matrix/XJ_clinic_raw')
combined_roc$sample <- rownames(train_topbins)
CAnonCA_testroc$sample <- rownames(test_topbins)
rownames(fragCAnontrain_topbins) %in%XJ_clinic_raw$sample

combined_roc <- left_join(combined_roc,XJ_clinic_raw,by='sample') 
CAnonCA_testroc <- left_join(CAnonCA_testroc,XJ_clinic_raw,by='sample')
### only cancer samples have the stage inform
combined_roc2 <- combined_roc[which(combined_roc$obs_num == 1),] # 84/88
combined_roc2[is.na(combined_roc2$stage),'stage_less'] = 'non'
CAnonCA_testroc2 <- CAnonCA_testroc[which(CAnonCA_testroc$obs_num == 1),] # 25/29
CAnonCA_testroc2[is.na(CAnonCA_testroc2$stage),'stage_less'] = 'non'

##### traing result
stage_training <- data.frame(Group = rep(c("I",'II','III','IV','NA'),each = 4),
                             Category = rep(c('cnv','tss','frag','combined'),5),
                             value = 1)
nrow(combined_roc2[which(combined_roc2$stage_less == 'I'),]) #12
nrow(combined_roc2[which(combined_roc2$stage_less == 'II'),]) #27
nrow(combined_roc2[which(combined_roc2$stage_less == 'III'),]) #27
nrow(combined_roc2[which(combined_roc2$stage_less == 'IV'),]) #6
nrow(combined_roc2[which(combined_roc2$stage_less == 'non'),]) #68


stage_res <- function(matrix,stage_list,model_num){
  combined_roc2_I <- matrix[which(matrix$stage_less %in% stage_list),]
  nsmall =  nrow(combined_roc2_I[which(combined_roc2_I[,model_num] > 0.4999999),])
  nall = nrow(combined_roc2_I)
  return(nall)
}

stage_res_posn <- function(matrix,stage_list,model_num){
  #matrix = combined_roc2
  #stage_list = c('I')
  #model_num = 5
  combined_roc2_I <- matrix[which(matrix$stage_less %in% stage_list),]
  nsmall =  nrow(combined_roc2_I[which(combined_roc2_I[,model_num] > 0.4999999),])
  nall = nrow(combined_roc2_I)
  return(nsmall/nall)
}

model_num_list_triain <- c(3,5,7,9)
for(i in 1:4){
  model_num  = model_num_list_triain[i]
  stage_training[i,3] = stage_res(combined_roc2,c('0','I'),model_num)
  stage_training[(i+4),3] = stage_res(combined_roc2,c('II'),model_num)
  stage_training[(i+8),3] = stage_res(combined_roc2,c('III'),model_num)
  stage_training[(i+12),3] = stage_res(combined_roc2,c('IV'),model_num)
  stage_training[(i+16),3] = stage_res(combined_roc2,c('non'),model_num)
}

for(i in 1:4){
  model_num  = model_num_list_triain[i]
  stage_training[i,4] = stage_res_posn(combined_roc2,c('0','I'),model_num)
  stage_training[(i+4),4] = stage_res_posn(combined_roc2,c('II'),model_num)
  stage_training[(i+8),4] = stage_res_posn(combined_roc2,c('III'),model_num)
  stage_training[(i+12),4] = stage_res_posn(combined_roc2,c('IV'),model_num)
  stage_training[(i+16),4] = stage_res_posn(combined_roc2,c('non'),model_num)
}

#### accuracy
stage_training$accuracy <- stage_training$V4/stage_training$value

###  stage plot ### training 
stage_training$Category <- factor(stage_training$Category,levels=c('tss','cnv','frag','combined'))

# 创建绘图对象
p <- ggplot(stage_training, aes(x = Category, y = V4, fill = Category))




p + geom_bar(stat = "identity", position = "dodge") +
  facet_grid(. ~ Group, space = "free") +
  labs(x = " ", y = "Accuracy", title = " ",col = 'Models (AUC)') +
  scale_fill_manual(values = c("combined" = '#c82423', "tss" = "#9ac9db", "frag" = "#f8ac8c",'cnv'="#2878b5"))+
  ggplot2::theme(
    panel.grid.major = element_line(colour = NA),
    panel.background = element_rect(fill = "transparent", colour = 'black'),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank()) + ggtitle(" ") +
  ggplot2::theme_bw() 

########################## test 
########################## testing CAnonCA_testroc2 from above
##########################

CAnonCA_testroc2 <- CAnonCA_testroc[which(CAnonCA_testroc$obs_num == 1),] # 25/29
CAnonCA_testroc2[is.na(CAnonCA_testroc2$stage),'stage_less'] = 'non'

##### traing result
stage_training <- data.frame(Group = rep(c("I",'II','III','IV','NA'),each = 4),
                             Category = rep(c('cnv','tss','frag','combined'),5),
                             value = 1)
nrow(combined_roc2[which(combined_roc2$stage_less == 'I'),]) #12
nrow(combined_roc2[which(combined_roc2$stage_less == 'II'),]) #27
nrow(combined_roc2[which(combined_roc2$stage_less == 'III'),]) #27
nrow(combined_roc2[which(combined_roc2$stage_less == 'IV'),]) #6
nrow(combined_roc2[which(combined_roc2$stage_less == 'non'),]) #68


stage_res <- function(matrix,stage_list,model_num){
  combined_roc2_I <- matrix[which(matrix$stage_less %in% stage_list),]
  nsmall =  nrow(combined_roc2_I[which(combined_roc2_I[,model_num] > 0.4999999),])
  nall = nrow(combined_roc2_I)
  return(nall)
}

stage_res_posn <- function(matrix,stage_list,model_num){

  combined_roc2_I <- matrix[which(matrix$stage_less %in% stage_list),]
  nsmall =  nrow(combined_roc2_I[which(combined_roc2_I[,model_num] > 0.4999999),])
  nall = nrow(combined_roc2_I)
  return(nsmall/nall)
}
model_num_list_triain <- c(4:7)

for(i in 1:4){
  model_num  = model_num_list_triain[i]
  stage_training[i,3] = stage_res(CAnonCA_testroc2,c('0','I'),model_num)
  stage_training[(i+4),3] = stage_res(CAnonCA_testroc2,c('II'),model_num)
  stage_training[(i+8),3] = stage_res(CAnonCA_testroc2,c('III'),model_num)
  stage_training[(i+12),3] = stage_res(CAnonCA_testroc2,c('IV'),model_num)
  stage_training[(i+16),3] = stage_res(CAnonCA_testroc2,c('non'),model_num)
}

for(i in 1:4){
  model_num  = model_num_list_triain[i]
  stage_training[i,4] = stage_res_posn(CAnonCA_testroc2,c('0','I'),model_num)
  stage_training[(i+4),4] = stage_res_posn(CAnonCA_testroc2,c('II'),model_num)
  stage_training[(i+8),4] = stage_res_posn(CAnonCA_testroc2,c('III'),model_num)
  stage_training[(i+12),4] = stage_res_posn(CAnonCA_testroc2,c('IV'),model_num)
  stage_training[(i+16),4] = stage_res_posn(CAnonCA_testroc2,c('non'),model_num)
}

###  stage plot ### training 
stage_training$Category <- factor(stage_training$Category,levels=c('tss','cnv','frag','combined'))

# 创建绘图对象
p <- ggplot(stage_training, aes(x = Category, y = V4, fill = Category))




p + geom_bar(stat = "identity", position = "dodge") +
  facet_grid(. ~ Group, space = "free") +
  labs(x = " ", y = "Accuracy", title = " ",col = 'Models (AUC)') +
  scale_fill_manual(values = c("combined" = '#c82423', "tss" = "#9ac9db", "frag" = "#f8ac8c",'cnv'="#2878b5"))+
  ggplot2::theme(
    panel.grid.major = element_line(colour = NA),
    panel.background = element_rect(fill = "transparent", colour = 'black'),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank()) + ggtitle(" ") +
  ggplot2::theme_bw() 

#

##### combine train & test

train_stage <- (combined_roc2[,c(8,3,5,7,9:ncol(combined_roc2))])

test_stage <- (CAnonCA_testroc2[,c(3:ncol(CAnonCA_testroc2))])
colnames(train_stage) <- colnames(test_stage)
colnames(train_stage)

stage_accu <- rbind(train_stage,test_stage)

stage_training <- data.frame(Group = rep(c("I",'II','III','IV','NA'),each = 4),
                             Category = rep(c('cnv','tss','frag','combined'),5),
                             value = 1)


model_num_list_triain <- c(2:5)

for(i in 1:4){
  model_num  = model_num_list_triain[i]
  stage_training[i,3] = stage_res(stage_accu,c('0','I'),model_num)
  stage_training[(i+4),3] = stage_res(stage_accu,c('II'),model_num)
  stage_training[(i+8),3] = stage_res(stage_accu,c('III'),model_num)
  stage_training[(i+12),3] = stage_res(stage_accu,c('IV'),model_num)
  stage_training[(i+16),3] = stage_res(stage_accu,c('non'),model_num)
}

for(i in 1:4){
  model_num  = model_num_list_triain[i]
  stage_training[i,4] = stage_res_posn(stage_accu,c('0','I'),model_num)
  stage_training[(i+4),4] = stage_res_posn(stage_accu,c('II'),model_num)
  stage_training[(i+8),4] = stage_res_posn(stage_accu,c('III'),model_num)
  stage_training[(i+12),4] = stage_res_posn(stage_accu,c('IV'),model_num)
  stage_training[(i+16),4] = stage_res_posn(stage_accu,c('non'),model_num)
}

###  stage plot ### training 
stage_training$Category <- factor(stage_training$Category,levels=c('tss','cnv','frag','combined'))

# 创建绘图对象
p <- ggplot(stage_training, aes(x = Category, y = V4, fill = Category))




p + geom_bar(stat = "identity", position = "dodge") +
  facet_grid(. ~ Group, space = "free") +
  labs(x = " ", y = "Accuracy", title = " ",col = 'Models (AUC)') +
  scale_fill_manual(values = c("combined" = '#c82423', "tss" = "#9ac9db", "frag" = "#f8ac8c",'cnv'="#2878b5"))+
  ggplot2::theme(
    panel.grid.major = element_line(colour = NA),
    panel.background = element_rect(fill = "transparent", colour = 'black'),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank()) + ggtitle(" ") +
  ggplot2::theme_bw() 

#






######################## validation result
####################################################################################
############################         GC         ##########################
####################################################################################
################ training
########

### only cancer samples have the source inform
combined_roc2 <- combined_roc[which(combined_roc$obs_num == 1),] # 84/88
CAnonCA_testroc2 <- CAnonCA_testroc[which(CAnonCA_testroc$obs_num == 1),] # 25/29

##### traing result
stage_training <- data.frame(Group = rep(c('Gastric Cancer','Colorectal cancer'),each = 4),
                             Category = rep(c('cnv','tss','frag','combined'),2),
                             value = 1)

stage_res <- function(matrix,stage_list,model_num){
  #matrix <- combined_roc2
  #stage_list = c('W')
  #model_num = 9
  
  combined_roc2_I <- matrix[which(matrix$source %in% stage_list),]
  nsmall =  nrow(combined_roc2_I[which(combined_roc2_I[,model_num] > 0.4999999),])
  nall = nrow(combined_roc2_I)
 
  return( nall)
}
stage_res_posn <- function(matrix,stage_list,model_num){
  #matrix <- combined_roc2
  #stage_list = c('W')
  #model_num = 9
  
  combined_roc2_I <- matrix[which(matrix$source %in% stage_list),]
  nsmall =  nrow(combined_roc2_I[which(combined_roc2_I[,model_num] > 0.4999999),])
  nall = nrow(combined_roc2_I)
  
  return( nsmall/nall)
}

 
model_num_list_triain <- c(3,5,7,9)
for(i in 1:4){
  model_num  = model_num_list_triain[i]
  stage_training[i,3] = stage_res(combined_roc2,c('W'),model_num)
  stage_training[i+4,3] = stage_res(combined_roc2,c('C'),model_num)
}

for(i in 1:4){
  model_num  = model_num_list_triain[i]
  stage_training[i,4] = stage_res_posn(combined_roc2,c('W'),model_num)
  stage_training[i+4,4] = stage_res_posn(combined_roc2,c('C'),model_num)
}

### accuracy
stage_training$accuracy <- stage_training$V4/stage_training$value

###  stage plot ### training 
stage_training$Category <- factor(stage_training$Category,levels=c('tss','cnv','frag','combined'))

# 创建绘图对象
p <- ggplot(stage_training, aes(x = Category, y = V4, fill = Category))




p + geom_bar(stat = "identity", position = "dodge") +
  facet_grid(. ~ Group, space = "free") +
  labs(x = " ", y = "Accuracy", title = " ",col = 'Models (AUC)') +
  scale_fill_manual(values = c("combined" = '#c82423', "tss" = "#9ac9db", "frag" = "#f8ac8c",'cnv'="#2878b5"))+
  ggplot2::theme(
    panel.grid.major = element_line(colour = NA),
    panel.background = element_rect(fill = "transparent", colour = 'black'),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank()) + ggtitle(" ") +
  ggplot2::theme_bw() 

########################## test 
##### traing result
stage_training <- data.frame(Group = rep(c('Gastric Cancer','Colorectal cancer'),each = 4),
                             Category = rep(c('cnv','tss','frag','combined'),2),
                             value = 1)

model_num_list_triain <- c(4:7)
for(i in 1:4){
  model_num  = model_num_list_triain[i]
  stage_training[i,3] = stage_res(CAnonCA_testroc2,c('W'),model_num)
  stage_training[i+4,3] = stage_res(CAnonCA_testroc2,c('C'),model_num)
}

for(i in 1:4){
  model_num  = model_num_list_triain[i]
  stage_training[i,4] = stage_res_posn(CAnonCA_testroc2,c('W'),model_num)
  stage_training[i+4,4] = stage_res_posn(CAnonCA_testroc2,c('C'),model_num)
}

### accuracy
stage_training$accuracy <- stage_training$V4/stage_training$value

###  stage plot ### training 
stage_training$Category <- factor(stage_training$Category,levels=c('tss','cnv','frag','combined'))

# 创建绘图对象
p <- ggplot(stage_training, aes(x = Category, y = V4, fill = Category))




p + geom_bar(stat = "identity", position = "dodge") +
  facet_grid(. ~ Group, space = "free") +
  labs(x = " ", y = "Accuracy", title = " ",col = 'Models (AUC)') +
  scale_fill_manual(values = c("combined" = '#c82423', "tss" = "#9ac9db", "frag" = "#f8ac8c",'cnv'="#2878b5"))+
  ggplot2::theme(
    panel.grid.major = element_line(colour = NA),
    panel.background = element_rect(fill = "transparent", colour = 'black'),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank()) + ggtitle(" ") +
  ggplot2::theme_bw() 


##### combine train & test

train_stage <- (combined_roc2[,c(8,3,5,7,9:ncol(combined_roc2))])

test_stage <- (CAnonCA_testroc2[,c(3:ncol(CAnonCA_testroc2))])
colnames(train_stage) <- colnames(test_stage)
colnames(train_stage)

stage_accu <- rbind(train_stage,test_stage)
stage_training <- data.frame(Group = rep(c('Gastric Cancer','Colorectal cancer'),each = 4),
                             Category = rep(c('cnv','tss','frag','combined'),2),
                             value = 1)


model_num_list_triain <- c(2:5)


for(i in 1:4){
  model_num  = model_num_list_triain[i]
  stage_training[i,3] = stage_res(stage_accu,c('W'),model_num)
  stage_training[i+4,3] = stage_res(stage_accu,c('C'),model_num)
}

for(i in 1:4){
  model_num  = model_num_list_triain[i]
  stage_training[i,4] = stage_res_posn(stage_accu,c('W'),model_num)
  stage_training[i+4,4] = stage_res_posn(stage_accu,c('C'),model_num)
}


###  stage plot ### training 
stage_training$Category <- factor(stage_training$Category,levels=c('tss','cnv','frag','combined'))

# 创建绘图对象
p <- ggplot(stage_training, aes(x = Category, y = V4, fill = Category))




p + geom_bar(stat = "identity", position = "dodge") +
  facet_grid(. ~ Group, space = "free") +
  labs(x = " ", y = "Accuracy", title = " ",col = 'Models (AUC)') +
  scale_fill_manual(values = c("combined" = '#c82423', "tss" = "#9ac9db", "frag" = "#f8ac8c",'cnv'="#2878b5"))+
  ggplot2::theme(
    panel.grid.major = element_line(colour = NA),
    panel.background = element_rect(fill = "transparent", colour = 'black'),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank()) + ggtitle(" ") +
  ggplot2::theme_bw() 

#

###############################################################
#######################fig 6b confusion table
draw_confusion_matrix <- function(cm) {
  
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="accuracy = 0.917(22/24)", ylab="", xaxt='n', yaxt='n')
  #title('CONFUSION MATRIX', cex.main=2)
  
  # create the matrix 
  rect(150, 430, 240, 370, col='navy')
  text(195, 435, 'Colorectal', cex=1.2)
  rect(250, 430, 340, 370, col='white')
  text(295, 435, 'Gastric', cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col='white')
  rect(250, 305, 340, 365, col='navy')# #F7AD50
  text(140, 400, 'Colorectal', cex=1.2, srt=90)
  text(140, 335, 'Gastric', cex=1.2, srt=90)
  
  # add in the cm results 
  
  text(195, 400, 0, cex=1.6, font=2, col='white')
  text(195, 335, 0, cex=1.6, font=2, col='black')
  text(295, 400, 2, cex=1.6, font=2, col='black')
  text(295, 335, 22, cex=1.6, font=2, col='white')
  res <- as.numeric(cm$table)
  text(195, 400, res[1], cex=1.6, font=2, col='black')
  text(195, 335, res[2], cex=1.6, font=2, col='black')
  text(295, 400, res[3], cex=1.6, font=2, col='black')
  text(295, 335, res[4], cex=1.6, font=2, col='black')
} 
### vali

tss_CAnon_trare <- data.frame(modeltss_WC$pred)

tss_CAnon_trare <- tss_CAnon_trare[which(tss_CAnon_trare$mtry == 1),]

tss_roc_result <- roc_mean_calculate(tss_CAnon_trare)
tss_rocplot <- roc_mean_calculate_matrix(tss_CAnon_trare)
tss_rocplot$pred_log <- 'C'
tss_rocplot[which(tss_rocplot$W > 0.49 ),'pred_log']='W'
tss_rocplot$obs_log <- 'C'
tss_rocplot[which(tss_rocplot$obs_num == 1 ),'obs_log']='W'

tss_rocplot <- cbind(tss_rocplot,model_clinic[which(model_clinic$sample %in%rownames(train_topbins_tssWC)),])


tss_rocplot$right <- 1
tss_rocplot[which(tss_rocplot$pred_log != tss_rocplot$obs_log ),'right']=0

tss_rocplot2 <- tss_rocplot[,c('stage_less','right','W','source')] 

tss_rocplot3 <- tss_rocplot2[which(tss_rocplot2$stage_less %in% c('II') ),]
tss_rocplot3 <- tss_rocplot2[is.na(tss_rocplot2 ),]

pred_WC_stage<- data.frame(satge = c('I','II','III','IV','NA'), accuracy= c(0.615,0.926, 0.8148,0.667,0.956))
###  stage plot ### training 
#stage_training$Category <- factor(stage_training$Category,levels=c('tss','cnv','frag','combined'))

# 创建绘图对象
p <- ggplot(pred_WC_stage, aes(x = satge, y = accuracy, fill = satge))

p + geom_bar(stat = "identity", position = "dodge") +
  labs(x = " ", y = "Accuracy", title = " ",col = 'Models (AUC)') +
  scale_fill_manual(values = c("I" = '#c82423', "II" = "#9ac9db", "III" = "#f8ac8c",'IV'="#2878b5",'NA'='navy'))+
  ggplot2::theme(
    panel.grid.major = element_line(colour = NA),
    panel.background = element_rect(fill = "transparent", colour = 'black'),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank()) + ggtitle(" ") +
  ggplot2::theme_bw() 


### train
confusion_table <- confusionMatrix(data = factor(tss_rocplot$pred_log), reference = factor(tss_rocplot$obs_log))
print(confusion_table)
draw_confusion_matrix(confusion_table)

### test
roc_data_tssWC


predWC <- predict(modeltss_WC,newdata = test_topbins_tssWC,type = "prob")
predWC$sample = rownames(predWC)
predWC <- left_join(predWC,model_clinic,'sample')

predWC <- predWC[,c(1:3,9)]
#predWC <- cbind(predWC,test_topbins_tssWC)
#tss_CAnon_trare <- data.frame(predWC[,1])

#tss_CAnon_trare <- tss_CAnon_trare[which(tss_CAnon_trare$mtry == 1),]

predWC$pred_log <- 'C'
predWC[which(predWC$W > 0.49 ),'pred_log']='W'
predWC$obs_log <- 'C'
predWC[which(predWC$source == 'W' ),'obs_log']='W'

predWC$right <- 1
predWC[which(predWC$pred_log != predWC$obs_log ),'right']=0

predWC_2 <- predWC[,c('stage_less','right','W','source')] 

predWC_3 <- predWC_2[which(predWC_2$stage_less %in% c('III') ),]
predWC_3 <- predWC_2[is.na(predWC_2 ),]

pred_WC_stage<- data.frame(satge = c('I','II','III','IV','NA'), accuracy= c(0.625,0.4,0.166,0.67,0.86))
###  stage plot ### training 
#stage_training$Category <- factor(stage_training$Category,levels=c('tss','cnv','frag','combined'))

# 创建绘图对象
p <- ggplot(pred_WC_stage, aes(x = satge, y = accuracy, fill = satge))

p + geom_bar(stat = "identity", position = "dodge") +
  labs(x = " ", y = "Accuracy", title = " ",col = 'Models (AUC)') +
  scale_fill_manual(values = c("I" = '#c82423', "II" = "#9ac9db", "III" = "#f8ac8c",'IV'="#2878b5",'NA'='navy'))+
  ggplot2::theme(
    panel.grid.major = element_line(colour = NA),
    panel.background = element_rect(fill = "transparent", colour = 'black'),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank()) + ggtitle(" ") +
  ggplot2::theme_bw() 


### train
confusion_table <- confusionMatrix(data = factor(predWC$pred_log), reference = factor(predWC$obs_log))
print(confusion_table)
draw_confusion_matrix(confusion_table)

######################## validation result
########################## 
stage_validation <- data.frame(Group = rep(c("I",'II','III','IV','NA'),each = 4),
                               Category = rep(c('cnv','tss','frag','combined'),5),
                               value = 1)

model_num_list_triain <- c(4:7)
for(i in 1:4){
  model_num  = model_num_list_triain[i]
  stage_validation[i,3] = stage_res(CAnonCA_testroc2,c('0','I'),model_num)
  stage_validation[(i+4),3] = stage_res(CAnonCA_testroc2,c('II'),model_num)
  stage_validation[(i+8),3] = stage_res(CAnonCA_testroc2,c('III'),model_num)
  stage_validation[(i+12),3] = stage_res(CAnonCA_testroc2,c('IV'),model_num)
  stage_validation[(i+16),3] = stage_res(CAnonCA_testroc2,c('non'),model_num)
}
###  stage plot ### training 
stage_training$Category <- factor(stage_training$Category,levels=c('cnv','tss','frag','combined'))

# 创建绘图对象
p <- ggplot(stage_training, aes(x = Category, y = V4, fill = Category))

# 添加柱状图层
p + geom_bar(stat = "identity", position = "dodge") +
  facet_grid(. ~ Group, space = "free") +
  labs(x = " ", y = "Accuracy", title = " ",col = 'Models (AUC)') +
  scale_fill_manual(values = c("combined" = '#c82423', "tss" = "#9ac9db", "frag" = "#f8ac8c",'cnv'="#2878b5"))+
  ggplot2::theme(
    panel.grid.major = element_line(colour = NA),
    panel.background = element_rect(fill = "transparent", colour = 'black'),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank())

###  stage plot ### validation
stage_validation
stage_validation$Category <- factor(stage_validation$Category,levels=c('cnv','tss','frag','combined'))

# 创建绘图对象
p <- ggplot(stage_validation, aes(x = Category, y = value, fill = Category))

# 添加柱状图层
p + geom_bar(stat = "identity", position = "dodge") +
  facet_grid(. ~ Group, space = "free") +
  labs(x = " ", y = "Accuracy", title = " ",col = 'Models (AUC)') +
  scale_fill_manual(values = c("combined" = '#c82423', "tss" = "#9ac9db", "frag" = "#f8ac8c",'cnv'="#2878b5"))+
  ggplot2::theme(
    panel.grid.major = element_line(colour = NA),
    panel.background = element_rect(fill = "transparent", colour = 'black'),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank()) + ggtitle(" ") +
  ggplot2::theme_bw() 


+coord_fixed(ratio = 10000000)


+ 
  ggplot2::theme(axis.text = element_text( color = "black", size = 15),
                 text = element_text(size = 25),  # 调整文字大小
                 plot.title = element_text(size = 20),  # 调整标题文字大小
                 plot.margin = unit(c(2, 2, 2, 2), "cm"),  # 调整边距
                 panel.grid.major=ggplot2::element_line(colour=NA),panel.grid.minor = ggplot2::element_blank())
