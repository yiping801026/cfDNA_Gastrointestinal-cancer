args=commandArgs(T)

orig_frag <- read.delim(args[1], header=FALSE)
### remove chrX chrY & only include the fragment length between 100-220
orig_frag_norm <- data.frame(orig_frag[-which(orig_frag$V1 %in% c('chrX','chrY')),c(1:2,102:223)])
# short frag: 100-150bp
orig_frag_norm$shrot_sum <- apply(orig_frag_norm[,3:53],1,sum)
# long frag: 150-220bp
orig_frag_norm$long_sum <- apply(orig_frag_norm[,54:124],1,sum)
orig_frag_norm <- data.frame(orig_frag_norm[,c(1:2,125:126)])
orig_frag_norm <- data.frame(orig_frag_norm)

frag_chr <- read.csv(args[4], row.names=1)

frag_result<- data.frame(chr=rep(frag_chr$chr,frag_chr$chr_number),
                         num=rep(frag_chr$chr_number,frag_chr$chr_number),
                         short=rep(frag_chr$chr_number,frag_chr$chr_number),
                         long=rep(frag_chr$chr_number,frag_chr$chr_number),
                         ratio=rep(frag_chr$chr_number,frag_chr$chr_number)
                         )

for(j in 1:22){
  N_j=frag_chr[j,2]

  N_end=sum(frag_chr[1:j,2])
  N_start=N_end-N_j+1

  one_chr_frag=orig_frag_norm[which(orig_frag_norm$V1 == frag_chr[j,'chr']),]
  one_chr_frag
  Nregion=round(nrow(one_chr_frag)/500)
  out_short <- rep(1,Nregion)
  out_long <- rep(1,Nregion)

  for (i in 1:Nregion) {
    start = (i-1)*500+1
    if (i!= Nregion){
      end = i*500-1
      out_short[i] <- sum(one_chr_frag[start:end,3])
      out_long[i] <- sum(one_chr_frag[start:end,4])
      } else{
        end = nrow(one_chr_frag)
        out_short[i] <- sum(one_chr_frag[start:end,3])
        out_long[i] <- sum(one_chr_frag[start:end,4])} }

  out_sum <- out_short/out_long
  frag_result[N_start:N_end,3] = out_short
  frag_result[N_start:N_end,4] = out_long
  frag_result[N_start:N_end,5] = out_sum

  }

setwd(args[3])
write.table(frag_result,args[2])

