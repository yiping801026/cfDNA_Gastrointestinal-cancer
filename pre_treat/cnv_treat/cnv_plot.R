args=commandArgs(T)
library(ggplot2)
### cnv plot function for individuals
cnv_individual_test <- function(cnv_test_matrix,he_raw_out,sample_name){
  
  ### pre-treat the matrix (remove the position name & add colnames)
  cnv_test_matrix <- cnv_test_matrix[,-2]
  colnames(cnv_test_matrix)[1] <- 'chr'
  raw_save_patient <- cnv_test_matrix
  
  ### for patients
  pa_raw <- data.frame(raw_save_patient[,2:ncol(raw_save_patient)])
  
  ### remove rows including '0' from healthy
  pa_zero_pos <- which(pa_raw==0)
  pa_zero_row <- as.integer(unique(sort(pa_zero_pos%%nrow(pa_raw))) )
  he_row <- as.numeric(gsub('X', '', rownames(he_raw_out)))
  pa_raw <- pa_raw[he_row,]
  
  ### normalize: calculate sum of patient && percentage to replace every cell
  t_pa_raw <- data.frame(t(pa_raw))
  t_pa_raw$sum <- apply(t_pa_raw,1,sum)
  pasum_col <- ncol(t_pa_raw)
  t_pa_raw_out <- t_pa_raw[,-pasum_col]/t_pa_raw[,pasum_col]
  pa_raw_out <- data.frame(t(t_pa_raw_out))
  
  ###  delete by median of healthy ones 
  per_pa_100CNV_wide <- data.frame(pa_raw_out/he_raw_out$median)
  
  per_pa_100CNV_wide$pos <- 1:nrow(per_pa_100CNV_wide)
  ################################################################
  ### log2-------healthy did log2 or not?
  ### without log2 range :  7.137759e-04 4.562386e+03
  ### with log2 range : -10.45224  12.15557   
  per_pa_100CNV_wide[,-ncol(per_pa_100CNV_wide)] <- log2(per_pa_100CNV_wide[,-ncol(per_pa_100CNV_wide)])
  write.csv(per_pa_100CNV_wide,paste0(sample_name,'per_pa_100CNV_wide.csv'))
  ################################################################
  ### wide to long
  per_pa_100CNV_wide$pos <- factor(per_pa_100CNV_wide$pos,levels = 1:nrow(per_pa_100CNV_wide))
  #summary(per_pa_100CNV_wide)
  N_samples_pa = ncol(per_pa_100CNV_wide)-1
  per_pa_100CNV_long <- tidyr::gather(per_pa_100CNV_wide,sample_ID,CNV_ratio,1:N_samples_pa,factor_key=TRUE)
  
  ### change position back to integer
  per_pa_100CNV_long$pos <- as.integer(per_pa_100CNV_long$pos )
  ### only use pos && CNV_ratio
  pa_plot <- per_pa_100CNV_long[,c(1,3)]
  
  ### plot 2
  
  p_he <-
    ggplot2::ggplot(pa_plot, ggplot2::aes(pa_plot$pos, pa_plot$CNV_ratio) ) +
    #xlim(0,27633)+
    ggplot2::ylim(-1,1)+ggplot2::ylab('log2 CNV ratio')+ ggplot2::xlab('bin')+
    #scale_fill_continuous(type = "viridis") +
    scico::scale_fill_scico(palette = "lajolla")+
    ggplot2::stat_density_2d(ggplot2::aes(fill = ..level..), geom = "polygon"
                             ,colour="white")
  
  p5 <- p_he+ggplot2::ggtitle(" ") + ggplot2::theme_bw() + ggplot2::theme(panel.grid=ggplot2::element_blank())
  
  pdf(file=paste0(sample_name,'cnv.pdf'),bg = "white")
  print(p5)
  dev.off()
  
  
}


### main 

#cnv <- read.delim("/Users/pingyi/Desktop/renke/antu/CNV100K",header=FALSE)
cnv <- read.delim(args[1],header=FALSE)
#he_raw_out <- read.csv("~/Desktop/cfDNA_treat/cnv_treat/he_raw_out.csv", row.names=1)
he_raw_out <- read.csv(args[2], row.names=1)
#file_list <- read.delim('/Users/pingyi/Desktop/renke/antu/cnv_file.txt',header=FALSE)
file_list <- read.delim(args[3],header=FALSE)
#result_out <- '/Users/pingyi/Desktop/renke/antu/cnv_result/'
result_out <- args[4]

for (i in 1:nrow(file_list)){
  print(file_list[i,])  
  sample_name=paste0(result_out,'/',as.character(file_list[i,]))
  sample_cnv_id = i+2
  cnv_individual_test(cnv[,c(1:2,sample_cnv_id)],he_raw_out,sample_name)
}

cnv_individual_test(cnv,he_raw_out,paste0(args[5],'/together'))

# Rscript /Users/pingyi/Desktop/cfDNA_treat/cnv_treat/cnv_plot.R "/Users/pingyi/Desktop/renke/antu/CNV100K" "~/Desktop/cfDNA_treat/cnv_treat/he_raw_out.csv" /Users/pingyi/Desktop/renke/antu/cnv_file.txt '/Users/pingyi/Desktop/renke/antu/cnv_result/'


