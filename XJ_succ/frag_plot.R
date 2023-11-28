### This script is used to plot for fragmentation profiles
### import library
library(readxl)
library(dplyr)
library(factoextra) 
library(ggplot2)
library(pROC)
library(stats)
library(caret)
library(reshape2)  
library(ggpubr)
library(ggstatsplot)

### box plot for the samples
### import the dataset && the clinic information 
#model_clinic <- read_excel("/Users/pingyi/Desktop/XJ_succ/output/data/model_clinic2.xlsx")
matrix_frag <- read.csv("/Users/pingyi/Desktop/XJ_succ/output/data/df_ratio_norm_sample.csv", row.names=1)
df_ratio_raw<- read.csv('/Users/pingyi/Desktop/XJ_succ/output/XJ_succ/pre_treat/frag_treat/frag_matrix_formodelall.csv', row.names=1)
model_clinic_stage <- read.csv("/Users/pingyi/Desktop/XJ_succ/output/data/model_clinic.csv",row.names =1)

df_ratio_raw <- df_ratio_raw[-c(184,146,115,340),]
matrix_frag_model=df_ratio_raw[,model_clinic_stage$sample_name]
dim(matrix_frag_model) # 558 450

### wide_to_long
frag_table <- matrix_frag_model
frag_long_table <- melt(frag_table)
colnames(frag_long_table) <- c('sample_name','frag_norm')
### frag_long_table boxplot for the Healthy vs Gastritis vs Cancer
#table(model_clinic$disease)
frag_long_table <- left_join(frag_long_table,model_clinic_stage,by = 'sample_name')

frag_long_table$type_less_subs <- frag_long_table$disease
frag_long_table[which(frag_long_table$type_less_subs == 'T'),'type_less_subs'] = 'Cancer'
frag_long_table[which(frag_long_table$type_less_subs == 'I'),'type_less_subs'] = 'Gastritis'
frag_long_table[which(frag_long_table$type_less_subs == 'H'),'type_less_subs'] = 'Healthy'


### add stages in the bottom 
unique(frag_long_table$stage_less)
frag_long_table_stage <- frag_long_table[!is.na(frag_long_table$stage_less),]
frag_long_table_stage$type_less_subs  = frag_long_table_stage$stage_less
unique(frag_long_table_stage$type_less_subs)
frag_long_table2 <- rbind(frag_long_table,frag_long_table_stage)
unique(frag_long_table2$type_less_subs)
### add locations in the bottom 
unique(frag_long_table$source)
frag_long_table_local <- frag_long_table[which(frag_long_table$disease == 'T'),]
frag_long_table_local$type_less_subs  =frag_long_table_local$source
frag_long_table_local[which(frag_long_table_local$type_less_subs == 'C'),
                      'type_less_subs']=  'Colorectal Cancer'
frag_long_table_local[which(frag_long_table_local$type_less_subs == 'W'),
                      'type_less_subs']=  'Gastric Cancer'
frag_long_table2 <- rbind(frag_long_table2,frag_long_table_local)
unique(frag_long_table2$type_less_subs)
### define colors for every groups
less_cbPalette <- c( '#FFB90F','#A1A9D0',
                    '#F0988C','#B883D4','#96CCCB','#CFEAF1','#F6CAE5',"brown3","#999999","navy")
### calculate the compare p values
compare_result <- compare_means(frag_norm ~ type_less_subs, data = frag_long_table2)

### to add p-values to the boxplot based on the unique(frag_long_table$type_less_subs)
my_comparisons <- list(c("Cancer","Gastritis"),c("Gastric Cancer","Colorectal Cancer"))
#, c("Cancer", "Healthy"),
#                       c("Gastritis", "Healthy"),c("0" ,"Healthy"),
#                       c("I" ,"Healthy"),c("II" ,"Healthy"),
#                       c("III" ,"Healthy"),c("IV" ,"Healthy"),
#                       c("Gastric Cancer" ,"Healthy"),c("Colorectal Cancer" ,"Healthy"))

frag_long_table2$type_less_subs = factor(frag_long_table2$type_less_subs,
         levels = c(  "Gastric Cancer" , "Colorectal Cancer",
                      "IV","III","II", "I" , "0"  , 
                      "Cancer" ,"Gastritis","Healthy"  ))
##################### boxplots for frag ##################### 
p <- ggplot(data=frag_long_table2,aes(x=type_less_subs,
                                     y=frag_norm,
                                     color=type_less_subs))+
  #geom_jitter(alpha=0.2,
  #            position=position_jitterdodge(jitter.width = 0.35, 
  #                                          jitter.height = 0, 
  #                                          dodge.width = 0.8))+
  
  geom_boxplot(alpha=0.2,width=0.45,
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+
  geom_violin(alpha=0.2,width=0.9,
              position=position_dodge(width=0.8),
              size=0.75)+
  stat_compare_means(label = 'p.signif',ref.group ='Healthy',label.y=1.1) + 
  stat_compare_means(comparisons = my_comparisons,label.y=0.8)+
  scale_color_manual(values = less_cbPalette)+
  theme_classic() +
  theme(legend.position="none") + 
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(face = 1,  size = 12,color = "black"),
        axis.text.y = element_text(face = 1,  size = 12,color = "black"))+ # x,y axis colors + text
  theme(text = element_text(size=16)) +ylab("Fragmentation ratio") + xlab("")
p + coord_flip()

##################### heatmap for frag features after treatment for model ##################### 

### import datasets
frag_heatmap_matrix <- read.csv("~/Desktop/XJ_succ/input_data/frag/frag_heatmap_matrix.csv",row.names = 1)
model_clinic <- model_clinic_stage
dim(frag_heatmap_matrix) 
### add annotation information to the sides -- samples(rows)
annotation_row <- data.frame( group = model_clinic$disease)
row.names(annotation_row) <- model_clinic$sample_name
annotation_row$stage <- model_clinic$stage_less
### change annotation names
annotation_row[which(annotation_row$group == 'I'),'stage'] = 'Gastritis' 
annotation_row[which(annotation_row$group == 'H'),'stage'] = 'Healthy' 
annotation_row[annotation_row=='T'] = 'Cancer'
annotation_row[annotation_row=='H'] = 'Healthy'
annotation_row[which(!annotation_row$group %in% c('Cancer','Healthy')),'group'] ='Gastritis'
annotation_row[is.na(annotation_row)] = 'NA'
annotation_row[annotation_row=='WI'] = 'Gastritis'

ann_colors = list(
  group = c(Healthy = "navy",Gastritis='#999999',`Cancer` ="brown3"),
  stage = c(Healthy ='navy',Gastritis='#999999',`0`='#F6CAE5',`I`='#CFEAF1',
            `II`='#96CCCB',`III`='#B883D4',`IV` = '#F0988C',`NA`='violet')
  
)


#raw_plot_data <- as.matrix(z_score_less[,-c(123:190,292,293)])
raw_plot_data <- data.frame(frag_heatmap_matrix)
#heat_plot <- rbind(raw_plot_data,raw_plot_data[grepl('WTXJ0P0162',rownames(raw_plot_data)),]) # 101 & 1(stage 0)
heat_plot <- raw_plot_data[model_clinic[which(model_clinic$stage_less == '0'),'sample'],] #(263)
heat_plot <- rbind(heat_plot,raw_plot_data[model_clinic[which(model_clinic$stage_less == 'I'),'sample'],]) # :264:283
heat_plot <- rbind(heat_plot,raw_plot_data[model_clinic[which(model_clinic$stage_less == 'II'),'sample'],]) #:315
heat_plot <- rbind(heat_plot,raw_plot_data[model_clinic[which(model_clinic$stage_less == 'III'),'sample'],]) #:348
heat_plot <- rbind(heat_plot,raw_plot_data[model_clinic[which(model_clinic$stage_less == 'IV'),'sample'],]) #:357
na_samples <- rownames(annotation_row[which(annotation_row$stage == 'NA'),])
heat_plot <- rbind(heat_plot,raw_plot_data[na_samples,])#:450
heat_plot <- rbind(heat_plot,raw_plot_data[model_clinic[which(model_clinic$disease == 'I'),'sample_name'],]) #(170:262)
heat_plot <- rbind(heat_plot,raw_plot_data[model_clinic[which(model_clinic$disease == 'H'),'sample_name'],]) #(170:262)

#######

library(dplyr)
library(pheatmap)
range(heat_plot)
### center + scale
heat_plot
preProcValues <- preProcess(heat_plot, method = c( "center",'scale'))
heat_plotTransformed <- predict(preProcValues, heat_plot)

p1 <- pheatmap(heat_plotTransformed,
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
               breaks = -4:4,
               color = colorRampPalette(c("darkblue", "white", "tomato2"))(9),
               #annotation_col=annotation_col
               
               #gaps_col=col_gap,
               gaps_row = c(188, 281),
               border_color = 'black'
               #border = TRUE
)

p1





