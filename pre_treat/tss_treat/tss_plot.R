library(caret)
library(readxl)
library(pheatmap)
args=commandArgs(T)
### import healthy ones
tss_healthy_raw <- read.csv(args[2], row.names=1)
tss_sample_raw <- read_excel(args[1])
#tss_sample_raw <- read_excel('/Users/pingyi/Desktop/renke/hos/BJRX/data/tssmatrix59.xlsx')
### combine samples && healthy together
tss_table <- cbind(tss_sample_raw,tss_healthy_raw)
tss_table_raw <- tss_table
############################# treat the tss matrix ############################# 
### remove zero && NAs
### change -1 into 1
tss_table_raw <- tss_table_raw[,3:ncol(tss_table_raw)]
tss_table_raw[tss_table_raw == -1] = 0
#### remove all 0s
tssMatrix_rowSum<-rowSums(tss_table_raw,na.rm=T)
tssMatrix_colSum<- colSums(tss_table_raw,na.rm=T)
tssMatrix_rm0<-tss_table_raw[tssMatrix_rowSum!=0,]
tssMatrix_rm0<-tssMatrix_rm0[,tssMatrix_colSum!=0]
dim(tssMatrix_rm0) # 92208   164
#### remove rows or cloumns with many NA and 0s
zero_matrix<-tssMatrix_rm0==0
zero_colSum<-colSums(zero_matrix,na.rm=T)
zero_rowSum<-rowSums(zero_matrix,na.rm=T)
na_matrix<-is.na(tssMatrix_rm0)
tss_naSum<-rowSums(na_matrix)
bin_naSum<-colSums(na_matrix)
tss_sampleNum<-dim(tssMatrix_rm0)[2]
tss_binNum<-dim(tssMatrix_rm0)[1]
tssMatrix_rmNA<-t(tssMatrix_rm0[(tss_naSum+zero_rowSum)<(tss_sampleNum/2),])#
dim(tssMatrix_rmNA)  # 164 90803
### center && scale for the dataset
preProcValues <- preProcess(tssMatrix_rmNA, method = c( "center",'scale'))
tssMatrix_trans <- predict(preProcValues, tssMatrix_rmNA)
#################################### before is for every dataset ####################################
#################################### after is for the several groups#################################
### add groups for tssMatrix_trans
tssMatrix_trans <- data.frame(tssMatrix_trans)
rownames(tssMatrix_trans)
#tssMatrix_trans$small_group <- c(rep('OT',38),rep('ZT',47),rep('H',79))
Nsamples = ncol(tss_sample_raw)-2
tssMatrix_trans$big_group <- c(rep('CA', Nsamples),rep('H',79))
bigLabels <- tssMatrix_trans$big_group
#smalllabels <- tssMatrix_trans$small_group
tssMatrix_trans <- tssMatrix_trans[,-ncol(tssMatrix_trans)]
### p-value of groups
binNum<-dim(tssMatrix_trans)[2] 
tout<-rep(1,binNum)

CA<-bigLabels=='CA'
for(i in 1:binNum)
{
  tout[i]<-t.test(tssMatrix_trans[CA,i], tssMatrix_trans[!CA,i])$p.value
}

tout_withID<-data.frame(cbind(tout,colnames(tssMatrix_trans)))
tout_withID[,1] <- as.numeric(tout_withID[,1])
range(tout_withID[,1]) #  3.189701e-22 9.999531e-01
dim(tout_withID) #  90803     2
tout_withID <- tout_withID[which(tout_withID[,1] < 0.05 ),]
dim(tout_withID) # 34736     2
tout_sort<-tout_withID[order(tout_withID[,1]),]
### use the p-value of two groups to select  the features
topbinNum<-1000 
topbins<-tout_sort[1:topbinNum,2]

tss_topbins<-as.data.frame(tssMatrix_trans)[,topbins]
dim(tss_topbins) # 164 1000
tss_topbins_CAH <- tss_topbins
write.csv(tss_topbins_CAH,file = paste0(args[3],'/tss_topbins_CAH.csv'))

### heatmap 
dim(tss_topbins)#164 1000
annotation_row <- data.frame( group = c(rep("Cancer",Nsamples),rep("Healthy",79)))
row.names(annotation_row) <- rownames(tss_topbins)
ann_colors = list(
  group = c(Cancer= "brown3",Healthy="navy")
)


p1 <- pheatmap(tss_topbins,
         #legend_breaks = 0:,
         cluster_rows = F,
         cluster_cols = T,
         scale = 'column',
         #scale = 'row',
         #scale = 'none',
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         show_colnames = F,
         show_rownames = F,
         cell_width = 0.1 ,
         cell_height = 0.1,
         breaks = -4:4,
         
         color = colorRampPalette(c("darkblue", "wheat3", "tomato2"))(8),
         #color = colorRampPalette(c("darkblue", "white", "red4"))(100),
         #annotation_col=annotation_col
         
         #gaps_col=col_gap,
         gaps_row = c(Nsamples),
         border_color = 'black'
         #border = TRUE
)

pdf(file=paste0(args[3],'/tss.pdf'),width=8,height=4,bg = "white")
print(p1)
dev.off()
