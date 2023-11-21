load('/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/tssCA.Rdata')
load('/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/tssWC.Rdata')
load('/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/tssWIH.Rdata')


####################################################################################
############################        CA vs  nonCA          ##########################
###################################
### GO analysis
tss_table_raw <- read.csv("~/Desktop/XJ_new/data/cnv/tss_table_raw.csv", row.names=1)
tss_colnames <- data.frame(row_id=rownames(tss_table_raw),genes=tss_table_raw$gene)
tss_colnames$row_id <- paste0('X',tss_colnames$row_id)
### find the gene names
tss_CAnon_features <- data.frame(row_id = rownames(modeltss_CAnonCA$finalModel$importance),
                                     importance=modeltss_CAnonCA$finalModel$importance)
tss_WIH_features <- data.frame(row_id = rownames(modeltss_WIH$finalModel$importance),
                                 importance=modeltss_WIH$finalModel$importance)
tss_WC_features <- data.frame(row_id = rownames(modeltss_WC$finalModel$importance),
                               importance=modeltss_CAnonCA$finalModel$importance)
tss_CAnon_features <- left_join(tss_CAnon_features,tss_colnames,by='row_id')
#tss_CAnon_features <- data.frame(row_id = rownames(tss_CAnonCA1000$finalModel$importance),
#                                 importance=tss_CAnonCA1000$finalModel$importance)

tss_WC_features <- left_join(tss_WC_features,tss_colnames,by='row_id')

tss_CAnonCA1000

library(clusterProfiler)
library(org.Hs.eg.db)
gene_list <- tss_CAnon_features$genes
write.csv(gene_list,'/Users/pingyi/Desktop/XJ_succ/output/gene_listCAnonCA.csv')
# 进行GO富集分析
go_result <- enrichGO(gene          = gene_list,
                      OrgDb         = org.Hs.eg.db, 
                      keyType       = "SYMBOL",       
                      ont           = "MF",           
                      pAdjustMethod = "BH",          
                      pvalueCutoff = 0.05)  
#CAnon249_MF <- go_result
CAnontss_MF <- data.frame(go_result@result)

#### order by the p-adjust
CAnontss_MF_top <- CAnontss_MF[order(CAnontss_MF$p.adjust),]
CAnontss_MF_top <- CAnontss_MF_top[1:10,]
CAnontss_MF_top$Generatio<- apply(data.frame(CAnontss_MF_top$GeneRatio),1,function(x)return(as.numeric(strsplit(x,'/')[[1]][1])/as.numeric(strsplit('2/118','/')[[1]][2])))
CAnontss_MF_top <- CAnontss_MF_top[order(-CAnontss_MF_top$Generatio),]
CAnontss_MF_top$Description)

CAnontss_MF_top$Description <- factor(CAnontss_MF_top$Description ,
                                      levels =c( "2-oxoglutarate-dependent dioxygenase activity"                                        
                                                 ,"calcium-dependent phospholipid binding"                                               
                                                 ,"steroid dehydrogenase activity"                                                       
                                                 ,"methyl-CpG binding"                                                                   
                                                 ,"ribosomal large subunit binding"                                                      
                                                 ,"oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor"
                                                 ,"calcium-dependent protein binding"                                                    
                                                 ,"cytokine activity"                                                                    
                                                 ,"signaling receptor activator activity"                                                
                                                 ,"receptor ligand activity"))
                                        

kegg_point=ggplot(CAnontss_MF_top,aes(x=Description,y=Generatio))+
  geom_point(aes(color=-log10(pvalue),size=Count),alpha=0.8)+
  scale_x_discrete(labels=function(x)stringr::str_wrap(x,width = 25))+
  coord_flip()+scale_size_continuous(range = c(3, 6)) +
  scale_color_gradient(low = "red",high = "purple")
kegg_point + ggplot2::theme(axis.text = element_text( color = "black", size = 10),
                            text = element_text(y='',size = 10),  # 调整文字大小
                            plot.title = element_text(size = 5),  # 调整标题文字大小
                            plot.margin = unit(c(0.1,0.1, 0.1, 0.1), "cm"))+coord_flip()

####################################################################################

############################        G vs  C             ##########################
###################################

  
### find the gene names
modeltss_WC_features <- data.frame(row_id = rownames(modeltss_WC$finalModel$importance),
                               importance=modeltss_WC$finalModel$importance)
modeltss_WC_features <- left_join(modeltss_WC_features,tss_colnames,by='row_id')
#tss_CAnon_features <- data.frame(row_id = rownames(tss_CAnonCA1000$finalModel$importance),
#                                 importance=tss_CAnonCA1000$finalModel$importance)

#tss_CAnon_features <- left_join(tss_CAnon_features,tss_colnames,by='row_id')

tss_CAnonCA1000

library(clusterProfiler)
library(org.Hs.eg.db)
gene_list <- modeltss_WC_features$genes
write.csv(gene_list,'/Users/pingyi/Desktop/XJ_succ/output/gene_listCAnonCA.csv')
# 进行GO富集分析
go_result <- enrichGO(gene          = gene_list,
                      OrgDb         = org.Hs.eg.db, 
                      keyType       = "SYMBOL",       
                      ont           = "CC",           
                      pAdjustMethod = "BH",          
                      pvalueCutoff = 0.05)  
#CAnon249_MF <- go_result
CAnontss_MF <- data.frame(go_result@result)

#### order by the p-adjust
CAnontss_MF_top <- CAnontss_MF[order(CAnontss_MF$p.adjust),]
CAnontss_MF_top <- CAnontss_MF_top[1:20,]
CAnontss_MF_top$Generatio<- apply(data.frame(CAnontss_MF_top$GeneRatio),1,function(x)return(as.numeric(strsplit(x,'/')[[1]][1])/as.numeric(strsplit('2/118','/')[[1]][2])))
CAnontss_MF_top <- CAnontss_MF_top[order(-CAnontss_MF_top$Generatio),]
CAnontss_MF_top$Description 




CAnontss_MF_top$Description <- factor(CAnontss_MF_top$Description ,
                                      levels =c( "calcium activated cation channel activity",
                                                 "calcium-activated potassium channel activity",
                                                 "beta-catenin binding", 
                                                 "ion gated channel activity",  
                                                 "heparin binding", 
                                                 "cytokine activity",
                                                 "actin filament binding",
                                                 "phosphoric ester hydrolase activity",
                                                 "signaling receptor activator activity", 
                                                 "receptor ligand activity")) 
                                      
                                      levels =c( "receptor ligand activity"   ,     
                                                 "signaling receptor activator activity",       
                                                 "phosphoric ester hydrolase activity", 
                                                 "actin filament binding",
                                                 "cytokine activity",
                                                 "heparin binding",                             
                                                 "ion gated channel activity",  
                                                 "beta-catenin binding",                        
                                                 "calcium-activated potassium channel activity",
                                                 "calcium activated cation channel activity"   ))


kegg_point=ggplot(CAnontss_MF_top,aes(x=Description,y=Generatio))+
  geom_point(aes(color=-log10(pvalue),size=Count),alpha=0.8)+
  scale_x_discrete(labels=function(x)stringr::str_wrap(x,width = 25))+
  coord_flip()+scale_size_continuous(range = c(3, 6)) +
  scale_color_gradient(low = "red",high = "purple")
kegg_point + ggplot2::theme(axis.text = element_text( color = "black", size = 10),
                            text = element_text(size = 10),  # 调整文字大小
                            plot.title = element_text(size = 5),  # 调整标题文字大小
                            plot.margin = unit(c(0.1,0.1, 0.1, 0.1), "cm"))+coord_flip()



tss_WC_features <- data.frame(row_id = rownames(modeltss_WC$finalModel$importance),
                              importance=modeltss_CAnonCA$finalModel$importance)




############################        WI vs  H            ##########################
###################################


### find the gene names
tss_WIH_features <- data.frame(row_id = rownames(modeltss_WIH$finalModel$importance),
                               importance=modeltss_WIH$finalModel$importance)
tss_WIH_features <- left_join(tss_WIH_features,tss_colnames,by='row_id')
#tss_CAnon_features <- data.frame(row_id = rownames(tss_CAnonCA1000$finalModel$importance),
#                                 importance=tss_CAnonCA1000$finalModel$importance)

#tss_CAnon_features <- left_join(tss_CAnon_features,tss_colnames,by='row_id')


gene_list <- tss_WIH_features$genes
#write.csv(gene_list,'/Users/pingyi/Desktop/XJ_succ/output/gene_listWIH.csv')
# 进行GO富集分析
go_result <- enrichGO(gene          = gene_list,
                      OrgDb         = org.Hs.eg.db, 
                      keyType       = "SYMBOL",       
                      ont           = "MF",           
                      pAdjustMethod = "BH",          
                      pvalueCutoff = 0.05)  
#CAnon249_MF <- go_result
CAnontss_MF <- data.frame(go_result@result)

#### order by the p-adjust
CAnontss_MF_top <- CAnontss_MF[order(CAnontss_MF$p.adjust),]
CAnontss_MF_top <- CAnontss_MF_top[1:10,]
CAnontss_MF_top$Generatio<- apply(data.frame(CAnontss_MF_top$GeneRatio),1,function(x)return(as.numeric(strsplit(x,'/')[[1]][1])/as.numeric(strsplit('2/118','/')[[1]][2])))
CAnontss_MF_top <- CAnontss_MF_top[order(CAnontss_MF_top$Generatio),]
CAnontss_MF_top$Description
CAnontss_MF_top$Description <- factor(CAnontss_MF_top$Description ,
                                      levels =c( "armadillo repeat domain binding",
                                                 "translation initiation factor binding",
                                                 "cullin family protein binding"  ,
                                                 "cAMP binding" ,
                                                 "protein N-terminus binding"  
                                                 ,"alcohol binding" ,
                                                 "GTPase binding" ,
                                                 "phosphatidylinositol binding"
                                                 ,"small GTPase binding","phospholipid binding"
                                                 ))
CAnontss_MF_top$Description <- factor(CAnontss_MF_top$Description ,
                                      levels =c( "phospholipid binding","small GTPase binding","phosphatidylinositol binding"         
                                                 ,"GTPase binding"                       
                                                 ,"alcohol binding"                      
                                                 ,"protein N-terminus binding"           
                                                 ,"cAMP binding"                         
                                                 ,"cullin family protein binding"        
                                                 ,"translation initiation factor binding"
                                                 ,"armadillo repeat domain binding"   ))


kegg_point=ggplot(CAnontss_MF_top,aes(x=as.factor(Description),y=Generatio))+
  geom_point(aes(color=-log10(pvalue),size=Count),alpha=0.8)+
  scale_x_discrete(labels=function(x)stringr::str_wrap(x,width = 25))+
  coord_flip()+scale_size_continuous(range = c(3, 6)) +
  scale_color_gradient(low = "red",high = "purple")
kegg_point + ggplot2::theme(axis.text = element_text( color = "black", size = 10),
                            text = element_text(size = 10),  # 调整文字大小
                            plot.title = element_text(size = 5),  # 调整标题文字大小
                            plot.margin = unit(c(0.1,0.1, 0.1, 0.1), "cm"))+coord_flip()



tss_WC_features <- data.frame(row_id = rownames(modeltss_WC$finalModel$importance),
                              importance=modeltss_CAnonCA$finalModel$importance)
#################################################################s
####################  heatmap of groups ########################
#####################################################################
#####################################################################
load(file='/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/tssCA.Rdata')
#alltss_zscore,tss_z_dhpbatchcut,DH_zscore_rmbatch_t,
#training_zscore_rmbatch,tss_z_pCAcut,train_topbins,test_topbins,
#modeltss_CAnonCA,pred,roc_data,
#frag_heatmap_matrix <- read.csv("~/Desktop/XJ_succ/input_data/frag/frag_heatmap_matrix.csv",row.names = 1)
model_clinic <- read.csv("~/Desktop/XJ_succ/input_data/clinic/model_clinic.csv",row.names =1)

dim(CAnon_matrix) # 450 301
CAnon_matrix<- rbind(train_topbins,test_topbins)

#### make heatmap

frag_heatmap_matrix <- CAnon_matrix[,-301]
dim() # 450 164
## sort
W_CA

raw_plot_data <- data.frame(frag_heatmap_matrix)
#heat_plot <- rbind(raw_plot_data,raw_plot_data[grepl('WTXJ0P0162',rownames(raw_plot_data)),]) # 101 & 1(stage 0)
heat_plot <- raw_plot_data[model_clinic[which(model_clinic$disease == 'T' & model_clinic$source == 'W'),'sample'],] #(263)

#heat_plot <- rbind(heat_plot,raw_plot_data[model_clinic[which(model_clinic$stage_less == '0'),'sample'],]) #(263)
heat_plot <- rbind(heat_plot,raw_plot_data[model_clinic[which(model_clinic$disease == 'T' & model_clinic$source == 'C'),'sample'],]) # :264:283
heat_plot <- rbind(heat_plot,raw_plot_data[model_clinic[which(model_clinic$disease == 'I'),'sample'],]) #:315
heat_plot <- rbind(heat_plot,raw_plot_data[model_clinic[which(model_clinic$disease == 'H'),'sample'],]) #:348

dim(heat_plot) # 450 300
model_clinic

#annotation_row <- data.frame( group = c(rep('Healthy',79),rep('Cancer',94)))
annotation_row <- data.frame( group = c(rep('Gastric Cancer',95),rep('Colorectal Cancer',93),rep('Inflammation',93),rep('Healthy',169)))

row.names(annotation_row) <- rownames(heat_plot)


ann_colors = list(
  group = c(`Gastric Cancer`="coral2" ,`Colorectal Cancer` = "brown3" ,Inflammation='salmon',`Healthy` ="navy")
  #stage = c(Healthy ='navy',`0 stage`='gold1',`I stage`='darkorange',`II stage`='chocolate',`III stage`='coral2',`IV stage` = 'brown4')
  
)


#colnames(CA_plot) <- apply(data.frame(names_CA),1,function(x) return(strsplit(x,'[.]')[[1]][1]))
### reorder by the stage 


#row.names(annotation_col) <- colnames(raw_plot_data)
library(pheatmap)
range(heat_plot)
### center + scale
heat_plot
preProcValues <- preProcess(heat_plot, method = c( "center",'scale'))
heat_plotTransformed <- predict(preProcValues, heat_plot)

p1 <- pheatmap(heat_plot,
               #legend_breaks = 0:,
               cluster_rows = F,
               cluster_cols = T,
               #scale = 'column',
               scale = 'row',
               #scale = 'none',
               annotation_row = annotation_row,
               annotation_colors = ann_colors,
               show_colnames = F,
               show_rownames = F,
               cell_width = 0.1 ,
               cell_height = 0.01,
               breaks = -4:5,
               color = colorRampPalette(c("darkblue", "white", "tomato2"))(8),
               
               #color = colorRampPalette(c("darkblue", "wheat3", "tomato2"))(8),
               #color = colorRampPalette(c("darkblue", "white", "red4"))(100),
               #annotation_col=annotation_col
               
               #gaps_col=col_gap,
               gaps_row = c(188, 281),
               border_color = 'black'
               #border = TRUE
)

p1


#################################################################s
####################  corr of C && W  ########################
#####################################################################
#####################################################################
load('/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/tssWC.Rdata')
load('/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/tssWIH.Rdata')
load('/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/tssCA.Rdata')


#save(train_topbins_tssWC,test_topbins_tssWC,
     modeltss_WC,pred_tssWC,roc_data_tssWC,
     file='/Users/pingyi/Desktop/XJ_succ/output/model_output_Rdata/tssWC.Rdata')



model_CW_matrix <- rbind(train_topbins,test_topbins)
model_CW_matrix <- model_CW_matrix[order(model_CW_matrix$group),]
model_CW_matrix_ls<- model_CW_matrix[,-301]

#dim(model_CW_matrix)#188 301
model_CW_matrix_lst <- data.frame(t(model_CW_matrix_ls))
#test_topbins_tssWC <- data.frame(t(test_topbins_tssWC[,-301]))

# 创建热图
heatmap(cor(model_CW_matrix_lst, method = "pearson"), 
        #Colv = NA, Rowv = NA,
        ColSideColors =  c(rep("brown3", 188), rep("orange", 262)),
        col = colorRampPalette(c("blue", "white", "brown3"))(50),  # 定义颜色范围
        symm = TRUE,  # 使矩阵对称
        main = "Pearson Correlation Heatmap",  # 图的标题
        xlab = " ", ylab = " "  # x和y轴标签
)


#### choose 15
model_CW_matrix_s <- cbind(t(train_topbins_tssWC),t(test_topbins_tssWC))
CW_matrix_cor <- cor(model_CW_matrix_s,method = "pearson")



#CAcor_matrix <- cor(frag_CA_norm, method = "pearson")
library(corrplot)
corrplot(CW_matrix_cor, type = "full",  tl.pos = "n",is.corr = TRUE)


ggcorr(CW_matrix_cor,type = 'full',nbreaks = 4, palette = "RdGy",tl.pos = "n"
       ,label = FALSE,geom=)

+ theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank()
)

##############################PCA of the C and W


PCA_plot_10PCs

model_CW_matrix <- rbind(train_topbins_tssWC,test_topbins_tssWC)
#model_CW_matrix <- model_CW_matrix[order(model_CW_matrix$group),]
model_CW_matrix_ls<- t(model_CW_matrix[,-301])
### import inhouse functions
source('/Users/pingyi/Desktop/XJ/XJ_paper/script/XJ_model/used_function.R')

PCA_plot_10PCs(model_CW_matrix_ls,model_CW_matrix$group,
               '/Users/pingyi/Desktop/XJ_succ/output/fig/fig3_tss/group_PCA/',20,20)

##############################PCA of the W and I

model_CW_matrix <- rbind(train_topbins_tssWIH,test_topbins_tssWIH)
#model_CW_matrix <- model_CW_matrix[order(model_CW_matrix$group),]
model_CW_matrix_ls<- t(model_CW_matrix[,-301])
### import inhouse functions
source('/Users/pingyi/Desktop/XJ/XJ_paper/script/XJ_model/used_function.R')

PCA_plot_10PCs(model_CW_matrix_ls,model_CW_matrix$group,
               '/Users/pingyi/Desktop/XJ_succ/output/fig/fig3_tss/group_PCAWIH/',150,150)

 
##############################################################
########### ggplot for C && W ###############
##############################################################
  library(reshape2)  
  library(tidyverse)
  library(hrbrthemes)
  library(viridis)
  library(forcats)
### CA
model_CW_matrix <- rbind(train_topbins_tssWC,test_topbins_tssWC)
model_CW_matrix <- rbind(train_topbins_tssWIH,test_topbins_tssWIH)

model_CW_matrix_ls <- melt(model_CW_matrix)

# for cnv plot fig2.2 about the variance differences between groups
library('ggstatsplot')
  p <- ggbetweenstats(
    data = model_CW_matrix_ls,
    x = group,
    y = value
  ) 
  
  p
  p+ ggplot2::scale_y_continuous(limits = c(0, 100))+
    ggplot2::theme(axis.text = element_text( color = "black", size = 10),
                   text = element_text(size = 15),  # 调整文字大小
                   plot.title = element_text(size = 10),  # 调整标题文字大小
                   plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")) 
  
#################### vol plot for CA nonCA###################
STADrna <- read.delim("~/Desktop/XJ_new/data/rna/STAD-胃癌.txt")
Crna <- read.delim("~/Desktop/XJ_new/data/rna/肠癌.txt")


library(ggplot2)
library(ggrepel)
############################# VOL plot function
VOL_plot <- function(result,med_name,result_data,control_name){
  ### make vol_data as input for plot making use padj
  #vol_data <- data.frame(gene=row.names(result), pval=-log10(result$padj), lfc=result$log2FoldChange)
  ### make vol_data as input for plot making use p-value
  vol_data <- data.frame(gene=STADrna$Gene.Symbol, pval=-log10(STADrna$adjp), lfc=STADrna$Log2.Fold.Change.)
  
  vol_data <- data.frame(gene=Crna$Gene.Symbol, pval=-log10(Crna$adjp), lfc=Crna$Log2.Fold.Change.)
  # remove zreo
  vol_data <- na.omit(vol_data)
  # set colors for up_regulated and down_regulated
  vol_data <- mutate(vol_data, color=case_when(
    vol_data$lfc >= 1 & vol_data$pval > -log10(0.05) ~ "UP",
    vol_data$lfc <= -1 & vol_data$pval > -log10(0.05) ~ "DOWN",
    vol_data$lfc > -1 & vol_data$lfc < 1 ~ "nonsignificant",
    vol_data$pval <= -log10(0.05) ~ "nonsignificant"))
  
  #### add the tssCAnon of 249 
  vol_data[which(vol_data$gene %in% tss_CAnon_features$genes),'color'] = 'tssCAnon'
  rnaCAtssgene <- data.frame(genes=vol_data[which(vol_data$color == 'tssCAnon'),'gene'] )
  #rnaCAtssgene <- left_join(rnaCAtssgene,tss_CAnon_249_features,by='genes')
  #### add the tssGC of 1000
  vol_data[which(vol_data$gene %in% tss_WC_features$genes),'color'] = 'tssCAnon'
  nrow(data.frame(genes=vol_data[which(vol_data$color == 'tssCAnon'),] ))
  
  vol_data[which(vol_data$gene %in% tss_WIH_features$genes),'color'] = 'tssCAnon'
  nrow(data.frame(genes=vol_data[which(vol_data$color == 'tssCAnon'),] ))
  
  
  
  #vol_data$color 
  #vol_data
  #tssCAnon
  
  # select genes we interested + top 40 genes
  vol_data$alpha = ifelse(vol_data$color == 'tssCAnon',1,0.8)
  vol_data$label = ifelse(vol_data$color == 'tssCAnon',as.character(vol_data$gene),'')
  nrow(vol_data[which(vol_data$color == 'tssCAnon'),])
  
  ### make VOL plot
  #vol_data <- vol_data_CAnontss_lab
  vol <- ggplot(vol_data, aes(x=lfc, y=pval,color=color ) )
  label <- data.frame(vol_data[which(vol_data$label != ''),])
  library(ggplot2)
  library(ggrepel)
  # change bands for VOL plot
  vol + 
    #ggtitle(label=paste(med_name,'_vs_',control_name," Volcano Plot",sep = ''), subtitle="Colored by fold-change direction") +
    geom_point(size=1, na.rm=T,aes(alpha =alpha)) +
    
    scale_color_manual(name="Regulated",
                       values=c(DOWN="#2878b5",UP="#c82423", tssCAnon='black',
                                nonsignificant="darkgray")) +
    theme_bw(base_size=14) +
    theme(legend.position="right") +
    xlab('log2 fold change') +
    ylab('-log10 adjust P value') +
    geom_hline(yintercept= -log10(0.05), colour="darkgrey",show.legend = TRUE) +
    geom_vline(xintercept=1, colour="darkgrey",show.legend = TRUE) +
    geom_vline(xintercept=-1, colour="darkgrey",show.legend = TRUE) +
    scale_y_continuous(limits = c(-1, 100))+
    scale_x_continuous(breaks = c(-5,-1,0,1,5),limits = c(-5, 5))+
    geom_text_repel(data =label , aes(x = lfc, y = pval, label = label),size = 4,max.overlaps =10000)
  #ggsave('_VOL.pdf',width = 20,height =10 )
  
  
  vol_data_CAnontss_lab <- vol_data
  #ggtitle(label=paste(med_name,'_vs_',control_name," Volcano Plot",sep = ''), subtitle="Colored by fold-change direction") +
  
  
}


save(tss_colnames,tss_CAnon_249_features,
     CAnon249_all, CAnon249_MF,tssGC_MF,tssGC_ALL,tssWIH_ALL,
     vol_data_GCtss_lab,vol_data_CAnontss_lab,
     file = "/Users/pingyi/Desktop/XJ_new/script/out_matrix/tss_plot.RData")


load("/Users/pingyi/Desktop/XJ_new/script/out_matrix/tss_plot.RData")

#######################################################################
################ PCA
tssMatrix_rmNA
colnames(tssMatrix_rmNA)
dim(tssMatrix_rmNA)
tssMatrix_rmNA2 <- data.frame(t(tssMatrix_rmNA))
tssMatrix_rmNA2$row_id <- paste0('X',rownames(tssMatrix_rmNA2))
### the rownames used in the CAnonCA (tss_WIH,tss_CAnonCA1000,tss_GC,tss_GC_features)
features <- rownames(data.frame(tss_WIH$finalModel$importance))
features <- rownames(data.frame(tss_CAnonCA1000$finalModel$importance))
features <- rownames(data.frame(tss_GC$finalModel$importance))

matrix <- tssMatrix_rmNA2[which(tssMatrix_rmNA2$row_id %in% features),]
dim(matrix)

library(factoextra) ### for PCA plot
library(ggplot2)
PCA_data_all <-  data.frame(t(matrix[,-264]))

dim(PCA_data_all)
#rownames(PCA_data_all)

#PCA_data_all <- PCA_data_all[c(1:60,129:291),]
#A_group_big <- c(rep('Cancer',117),rep('Inflammation',67),rep('Healthy',79))

A_group_big <- c(rep('Colorectal Cancer',55),rep('Gastric Cancer',62),rep('Inflammation',67),rep('Healthy',79))
PCA_data_all$A_group_big <- A_group_big

## WI
PCA_data_all <- PCA_data_all[118:263,]
## GC
PCA_data_all <- PCA_data_all[1:117,]

pca_F_result<- prcomp(PCA_data_all[,1:1000])

fviz_pca_ind(pca_F_result, label="none", habillage=PCA_data_all$A_group,
             addEllipses=T, ellipse.level=0.8,
             palette = c('brown3','#FFB90F','navy',"#999999"))+
  ggtitle(" ") +
  ggplot2::theme_bw() +
  #ggplot2::ylim(-20,30)+
  #ggplot2::xlim(-10,10)+
  ggplot2::theme(panel.grid.major=ggplot2::element_line(colour=NA),panel.grid.minor = ggplot2::element_blank())

###### heatmap 
dim(PCA_data_all)
annotation_row <- data.frame( group = c(rep('Cancer',117),rep('Inflammation',67),rep('Healthy',79)))
row.names(annotation_row) <- rownames(PCA_data_all)
#annotation_row$stage <- c(c(rep('Healthy',79),rep('0 stage',1),rep('I stage',22),rep('II stage',34),rep('III stage',29),
#rep('IV stage',8)))
### change the color of bar 

ann_colors = list(
  group = c(Healthy = "navy",Cancer="brown3",Inflammation='#999999')
)
library(pheatmap)
#PCA_data_all2 <- log2(PCA_data_all)


pheatmap(PCA_data_all,
         #legend_breaks = 0:,
         cluster_rows = T,
         cluster_cols = F,
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
         gaps_row = c(117),
         border_color = 'black'
         #border = TRUE
)


features <- rownames(data.frame(tss_GC$finalModel$importance))
matrix <- tssMatrix_rmNA2[which(tssMatrix_rmNA2$row_id %in% features),]
dim(matrix)
PCA_data_all <-  data.frame(t(matrix[,-264]))

dim(PCA_data_all)
PCA_data_all <- PCA_data_all[1:117,]
annotation_row <- data.frame( group = c(rep('Colorectal Cancer',55),rep('Gastric Cancer',62)))
row.names(annotation_row) <- rownames(PCA_data_all)
#annotation_row$stage <- c(c(rep('Healthy',79),rep('0 stage',1),rep('I stage',22),rep('II stage',34),rep('III stage',29),
#rep('IV stage',8)))
### change the color of bar 

ann_colors = list(
  group = c(`Colorectal Cancer` = "brown3",`Gastric Cancer`="#FFB90F")
)
library(pheatmap)

pheatmap(PCA_data_all,
         #legend_breaks = 0:,
         cluster_rows = T,
         cluster_cols = F,
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
         gaps_row = c(55),
         border_color = 'black'
         #border = TRUE
)




