#!/bin/bash

motu=$1
output_dir=$2
DnoisE_dir=$3
motus_tab_dir=$4
cores=$5
last_col=$(head -n 1 ${motus_tab_dir}/${motu} | grep -o "sample:" | wc -l)

if [ -f "${output_dir}/${motu}_Adcorr_denoised_ratio_d.csv" ]
then 
	echo "${output_dir}/${motu}_Adcorr_denoised_ratio_d.csv exists"
else
	python3 ${DnoisE_dir}/DnoisE.py --csv_input ${motus_tab_dir}/${motu} --csv_output ${output_dir}/${motu} -c ${cores} -s 4 -z ${last_col} -a 4 -n count -e 0.4881,0.2360,1.0078 -y 2>${output_dir}/${motu}_Adcorr_progress
fi