#' frag_ratio
#'
#' NO2.STEP
#' @param orig_file_path original file path for all frag.csv files

#'
#' @return df_ratio
#' @export
#' @importFrom utils read.csv
#' @examples
args=commandArgs(T)

frag_ratio <- function(orig_file_path){
  #orig_file_path = '/Users/pingyi/Desktop/XJ_new/frag_tmp/tmp'
  
  ### read in the csv files
  file_names_path <- dir(orig_file_path, pattern = "*.csv", recursive = F, full.names = T)

  ### DataFrame to store all ratio
  df_ratio_all <- read.csv(file_names_path[1],sep="")

  for (i in 2:length(file_names_path)) {
    df_ratio_all <- cbind(df_ratio_all, read.csv(file_names_path[i],sep=""))}

  ### change the file_names
  file_names <- dir(orig_file_path, pattern = "*.csv", recursive = F, full.names = F)

  file_names_less <- c(1:length(file_names))
  for (i in 1:length(file_names)){
    m <- strsplit(file_names[i],'[.]')[[1]][1]
    file_names_less[i] <- strsplit(m,'_')[[1]][1]}

  ### remove ratio is NA
  df_ratio_all <- df_ratio_all[!is.na(df_ratio_all$ratio),]

  ### get the ratios
  df_ratio <- df_ratio_all[,seq(from=5, to=ncol(df_ratio_all), by=5)]
  colnames(df_ratio) <- file_names_less

  ### without pos
  df_ratio <- data.frame(df_ratio)
  rownames(df_ratio) <- 1:nrow(df_ratio)

  ###add the chromosome infom
  df_ratio <- cbind(df_ratio,df_ratio_all$chr)

  ### scale for individuals
  N_samples = ncol(df_ratio)-1
  df_ratio[nrow(df_ratio)+1,1:N_samples] <-  apply(df_ratio[,1:N_samples],2,mean)

  ### subtract the mean of individual
  N_fragbins = nrow(df_ratio)-1
  for(i in 1:N_samples){
    df_ratio[1:N_fragbins,i] <-  (df_ratio[1:N_fragbins,i] - df_ratio[N_fragbins+1,i])
  }

  ### remove the mean
  df_ratio <- df_ratio[-nrow(df_ratio),]

  

  return(df_ratio)
}

### main
df_ratio = frag_ratio(args[1])
write.csv(df_ratio,args[2])
