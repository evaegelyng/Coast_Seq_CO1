#!/bin/bash

i=$1
MOTUS2RUN=$2
output_dir=$3

motu=$(sed "${i}q;d" ${MOTUS2RUN})
  		if [ $i == 1 ]
   		then
    			head -n 1 ${output_dir}${motu}_Adcorr_denoised_ratio_d.csv > ${output_dir}ESV_Adcorr.csv
    			sed -i 's/id/motu,id/g' ${output_dir}ESV_Adcorr.csv
   		fi
  		awk 'NR>1' ${output_dir}${motu}_Adcorr_denoised_ratio_d.csv | awk -v var="$motu," '{print var $0;}' >>${output_dir}ESV_Adcorr.csv
