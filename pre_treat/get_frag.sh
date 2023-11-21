### the original frag files : *chrBinFragLenDist, the tmp_dir used to store, file_list is the file names
patient_dir=$1
input_dir=$patient_dir/chr
output_dir=$patient_dir/results
tss_input=$patient_dir/tssmatrix.xlsx

### R script for frag calculate && plot && used files
R_script=/home/pingyi/script_pool/cfdna_treat/frag_treat/get_frag_ratio.R
script_data=/home/pingyi/script_pool/cfdna_treat/frag_treat/frag_chr.csv
R_script_norm=/home/pingyi/script_pool/cfdna_treat/frag_treat/frag_ratio.R
R_script_raw=/home/pingyi/script_pool/cfdna_treat/frag_treat/frag_ratio_raw.R
R_script_plot=/home/pingyi/script_pool/cfdna_treat/frag_treat/frag_plot.R
R_script_togetherplot=/home/pingyi/script_pool/cfdna_treat/frag_treat/frag_line_plot.R
frag_healthytable_raw=/home/pingyi/script_pool/cfdna_treat/frag_treat/frag_healthytable_raw.csv
frag_healthytable_norm=/home/pingyi/script_pool/cfdna_treat/frag_treat/frag_healthytable_norm.csv

### R script for tss plot && used files
R_script_tssplot=/home/pingyi/script_pool/cfdna_treat/tss_treat/tss_plot.R
tss_helathy=/home/pingyi/script_pool/cfdna_treat/tss_treat/tss_healthy_raw.csv

### output directories && files
tmp_dir=$output_dir/chr
plot_dir=$output_dir/chrplot
file_list=$input_dir/file.txt

frag_out_matrix=$output_dir/frag_out_matrix.csv
frag_out_matrix_raw=$output_dir/frag_out_matrix_raw.csv

### mkdir for results
mkdir $output_dir
cd $output_dir
mkdir $tmp_dir
mkdir $plot_dir
wait
#############################################################################
############### STEP1. pre_treat the fragmentation profile ##################
#############################################################################
cd $input_dir
ls *_chrBinFragLenDist > file.txt
wait 

### calculate the frag ratio for all samples;results store in tmp_dir.; 'script_data':number of bins of different chromosome
cat $file_list | while read file_name
do
    echo $file_name
    Rscript $R_script $input_dir/"$file_name" $tmp_dir/"${file_name}.csv" $tmp_dir $script_data
done
### norm by minus individual mean
wait
Rscript $R_script_raw $tmp_dir $frag_out_matrix_raw
Rscript $R_script_norm $tmp_dir $frag_out_matrix
#############################################################################
############### STEP2. plots for fragmentation profiles #####################
#############################################################################
### plot the frag line for individuals 
wait
Rscript $R_script_plot $frag_out_matrix $file_list $plot_dir
wait
### plot the frag line for all 
cd $output_dir
Rscript $R_script_togetherplot $frag_out_matrix $frag_out_matrix_raw $frag_healthytable_raw $frag_healthytable_norm

#############################################################################
############### STEP3. plots for tss (!not used in this paper) ##############
#############################################################################

################ plot for tss
Rscript $R_script_tssplot $tss_input $tss_helathy $output_dir












