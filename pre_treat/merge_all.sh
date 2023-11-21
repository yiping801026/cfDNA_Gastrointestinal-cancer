patient_dir=$1
input_dir=$patient_dir/cnv

### output directories && files
output_dir=$patient_dir/results
tmp_dir=$output_dir/cnv
file_list=$input_dir/file.txt

### R script for cnv calculate && plot && used files
one_input_file=/home/pingyi/script_pool/cfdna_treat/cnv_treat/cnv_firsttwo
R_script_plot=/home/pingyi/script_pool/cfdna_treat/cnv_treat/cnv_plot.R
he_raw_out=/home/pingyi/script_pool/cfdna_treat/cnv_treat/he_matrix_frag.csv

### make files for output
cd $input_dir
ls *CNV_100000 > file.txt
wait 

cd ..
mkdir ./results
mkdir ./results/cnv

### step 0, remove chrX +chrY for all CNV


cat $file_list | while read file_name
do
    echo $file_name
    cd $input_dir
    awk '{if($1!="chrY"&&$1!="chrX")print$3}' $file_name > $tmp_dir/$file_name
done

### step 1, first 2 columns of chr1:chr22 

cp $one_input_file $tmp_dir
wait

### step 2, merge together 
wait
cd $tmp_dir
paste cnv_firsttwo *100000 > "CNV100K"

### step 3, remove tmp files
rm *100000 
rm cnv_firsttwo
wait

### step 4, norm by using percentage as the input data & devide by the healthy ones & plot for cnv 
Rscript $R_script_plot $tmp_dir/CNV100K $he_raw_out $file_list $tmp_dir $output_dir



