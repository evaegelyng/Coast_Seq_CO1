#!/bin/bash

i=$1
MOTUS2RUN=$2
output_dir=$3
DnoisE_dir=$4
motus_tab_dir=$5
cores=$6

var=$(sed "${i}q;d" ${MOTUS2RUN})
  		if [ -f "${output_dir}/${var}_Adcorr_denoised_ratio_d.csv" ]
   		then 
    			echo "${output_dir}/${var}_Adcorr_denoised_ratio_d.csv exists" &
   		else
    			python3 ${DnoisE_dir}/DnoisE.py --csv_input ${motus_tab_dir}/${var} --csv_output ${output_dir}/${var} -c ${cores} -s 4 -z 78 -a 4 -n count -e 0.4902,0.2385,1.0088 -y 2>${output_dir}/${var}_Adcorr_progress &
     	fi