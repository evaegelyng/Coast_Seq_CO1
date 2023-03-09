#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 8G
#SBATCH -c 1
#SBATCH -t 04:00:00

input_dir="tmp/output_Ad_corr"
output_dir="tmp"
MOTUS2RUN="results/old_selected.txt"
lines=$(wc -l ${MOTUS2RUN} | cut -f1 -d ' ')

for i in $(seq 1 ${lines})
do
    motu=$(sed "${i}q;d" ${MOTUS2RUN})
    if [ $i == 1 ]
    then
	    head -n 1 ${input_dir}/${motu}_Adcorr_denoised_ratio_d.csv > ${output_dir}/ESV_Adcorr.csv
	    sed -i 's/id/motu,id/g' ${output_dir}/ESV_Adcorr.csv
    fi
    awk 'NR>1' ${input_dir}/${motu}_Adcorr_denoised_ratio_d.csv | awk -v var="$motu," '{print var $0;}' >>${output_dir}/ESV_Adcorr.csv
done