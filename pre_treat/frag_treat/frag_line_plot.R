library(ggplot2)
args=commandArgs(T)

frag_line_plot <- function(frag_matrix,N_health_samples,N_cancer_samples,type,sample_name){
  #frag_matrix = frag_out_matrix
  #N_health_samples = 1
  #N_cancer_samples = 59
  #type = 'non'

  if(type == 'control'){
    start_pos = 1
    end_pos = N_health_samples
    line_color = 'navy'
    plot_title = "Fragment ratio pattern of control samples"} else {
      start_pos = 1 + N_health_samples
      end_pos = N_health_samples + N_cancer_samples
      #line_color = 'brown3'
      line_color = 'goldenrod3'
      plot_title = "Fragment ratio pattern of test samples"}
  
  ### plot the samples
  df_ratio_plot2 <- data.frame(frag_matrix[,start_pos:end_pos])
  
  df_ratio_plot2$pos <- 1:nrow(df_ratio_plot2)
  df_ratio_plot2 <- data.frame(reshape2::melt(df_ratio_plot2,id="pos"))
  colnames(df_ratio_plot2)[2:3] <- c('sample','ratio')
  
  df_ratio_plot2$pos <- as.integer(df_ratio_plot2$pos)
  df_ratio_plot2$sample <- as.character(df_ratio_plot2$sample)
  
  df_ratio_plot2 <- df_ratio_plot2[which(! df_ratio_plot2$pos %in% c(184,146,115,340)),]
  
  p2 <- ggplot2::ggplot(data = df_ratio_plot2,ggplot2::aes(x=pos,y=ratio,group = sample,color=sample))+
    ggplot2::ggtitle(plot_title)+
    #ggplot2::xlab('Bin positions')+
    #ggplot2::ylab('Normalized fragment ratio')+
    ggplot2::ylim(0.0,0.8)+
    #ggplot2::ylim(-0.2,0.2)+
    #ggplot2::geom_line(color=line_color)+
    ggplot2::geom_line(color=line_color)+
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major=ggplot2::element_line(colour=NA),
                   panel.background = ggplot2::element_rect(fill = "transparent",colour = NA),
                   plot.background = ggplot2::element_rect(fill = "transparent",colour = NA),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.position = c(.075,.915),
                   legend.box.background = ggplot2::element_rect(color="black"))
  
  p3 <- p2 + scale_x_continuous(breaks = c(1,48,96,136,174,210,244,276,305,331,358,385,412,432,450,467,484,
                                           501,517,529,542,550,557.7 #23 start & end position of chromosome
                                           ,24.5,72.0, 116.0, 155.0 ,192.0, 227.0, 260.0 ,290.5, 318.0, 344.5,
                                           371.5, 398.5, 422.0 ,441.0 ,458.5 ,475.5, 492.5 ,509.0,
                                           523.0, 535.5, 546.0 ,554.0),    # for the x-axis names 
                                labels = c(rep(' ',23 ),1:22),expand=c(0,0))+ #rep(paste0("chr", 1:22), each = 1)
    theme(axis.line = element_line(colour = "black"),
          axis.text.x = element_text(face = 1,  size = 6.5,color = "black"),
          axis.text.y = element_text(face = 1,  size = 6.5,color = "black"))+ # x,y axis colors + text
    theme(text = element_text(size=10)) +ylab("Normalized fragment ratio") + xlab("chromosome")
  
  
  pdf(file=paste0(sample_name,'.pdf'),width=10,height=4,bg = "white")
  print(p3)
  dev.off()

}

frag_line_plot_norm <- function(frag_matrix,N_health_samples,N_cancer_samples,type,sample_name){
  #frag_matrix = frag_out_matrix
  #N_health_samples = 1
  #N_cancer_samples = 59
  #type = 'non'

  if(type == 'control'){
    start_pos = 1
    end_pos = N_health_samples
    line_color = 'navy'
    plot_title = "Fragment ratio pattern of control samples"} else {
      start_pos = 1 + N_health_samples
      end_pos = N_health_samples + N_cancer_samples
      #line_color = 'brown3'
      line_color = 'goldenrod3'
      plot_title = "Fragment ratio pattern of test samples"}
  
  ### plot the samples
  df_ratio_plot2 <- data.frame(frag_matrix[,start_pos:end_pos])
  
  df_ratio_plot2$pos <- 1:nrow(df_ratio_plot2)
  df_ratio_plot2 <- data.frame(reshape2::melt(df_ratio_plot2,id="pos"))
  colnames(df_ratio_plot2)[2:3] <- c('sample','ratio')
  
  df_ratio_plot2$pos <- as.integer(df_ratio_plot2$pos)
  df_ratio_plot2$sample <- as.character(df_ratio_plot2$sample)
  
  df_ratio_plot2 <- df_ratio_plot2[which(! df_ratio_plot2$pos %in% c(184,146,115,340)),]
  
  p2 <- ggplot2::ggplot(data = df_ratio_plot2,ggplot2::aes(x=pos,y=ratio,group = sample,color=sample))+
    ggplot2::ggtitle(plot_title)+
    #ggplot2::xlab('Bin positions')+
    #ggplot2::ylab('Normalized fragment ratio')+
    ggplot2::ylim(-0.2,0.2)+
    #ggplot2::ylim(-0.2,0.2)+
    #ggplot2::geom_line(color=line_color)+
    ggplot2::geom_line(color=line_color)+
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major=ggplot2::element_line(colour=NA),
                   panel.background = ggplot2::element_rect(fill = "transparent",colour = NA),
                   plot.background = ggplot2::element_rect(fill = "transparent",colour = NA),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.position = c(.075,.915),
                   legend.box.background = ggplot2::element_rect(color="black"))
  
  p3 <- p2 + scale_x_continuous(breaks = c(1,48,96,136,174,210,244,276,305,331,358,385,412,432,450,467,484,
                                           501,517,529,542,550,557.7 #23 start & end position of chromosome
                                           ,24.5,72.0, 116.0, 155.0 ,192.0, 227.0, 260.0 ,290.5, 318.0, 344.5,
                                           371.5, 398.5, 422.0 ,441.0 ,458.5 ,475.5, 492.5 ,509.0,
                                           523.0, 535.5, 546.0 ,554.0),    # for the x-axis names 
                                labels = c(rep(' ',23 ),1:22),expand=c(0,0))+ #rep(paste0("chr", 1:22), each = 1)
    theme(axis.line = element_line(colour = "black"),
          axis.text.x = element_text(face = 1,  size = 6.5,color = "black"),
          axis.text.y = element_text(face = 1,  size = 6.5,color = "black"))+ # x,y axis colors + text
    theme(text = element_text(size=10)) +ylab("Normalized fragment ratio") + xlab("chromosome")
  
  
  pdf(file=paste0(sample_name,'.pdf'),width=10,height=4,bg = "white")
  print(p3)
  dev.off()

}

frag_healthytable_raw <- read.csv(args[3], row.names=1)
frag_healthytable_norm <- read.csv(args[4], row.names=1)

df_ratio_norm <- read.csv(args[1],row.names=1)
df_ratio_raw <- read.csv(args[2],row.names=1)

frag_CAH_raw <- cbind(frag_healthytable_raw,df_ratio_raw[,-ncol(df_ratio_raw)])
frag_CAH_norm <- cbind(frag_healthytable_norm,df_ratio_norm[,-ncol(df_ratio_norm)])


frag_line_plot(frag_matrix = frag_CAH_raw, N_health_samples= ncol(frag_healthytable_raw),
               N_cancer_samples = ncol(df_ratio_raw)-1 ,type='non',sample_name = 'fragraw_line_plot')
frag_line_plot_norm(frag_matrix = frag_CAH_norm, N_health_samples= ncol(frag_healthytable_norm),
               N_cancer_samples = ncol(df_ratio_norm)-1 ,type='non',sample_name = 'fragnorm_line_plot')

frag_line_plot(frag_matrix = frag_CAH_raw, N_health_samples= ncol(frag_healthytable_raw),
               N_cancer_samples = ncol(df_ratio_raw)-1 ,type='control',sample_name = 'fragrawcontrol_line_plot')

frag_line_plot_norm(frag_matrix = frag_CAH_norm, N_health_samples= ncol(frag_healthytable_norm),
               N_cancer_samples = ncol(df_ratio_norm)-1 ,type='control',sample_name = 'fragnormcontrol_line_plot')
