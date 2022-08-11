#!/bin/bash

motus_tab_dir=$1
i=$2
MOTUS2RUN=$3
input_file=$4
motus_dir=$5

var=$(sed "${i}q;d" ${MOTUS2RUN})
sed "1q;d" ${input_file} >${motus_tab_dir}${var}
mkdir -p ${motus_tab_dir}/${var}_dir
temporal_dir=${motus_tab_dir}/${var}_dir/
#cd ${temporal_dir}
split -l200 -d ${motus_dir}/${var}
for z in ${temporal_dir}/*
do
        
	    linesmotu=$(sed -e ':a;N;$!ba;s/\n/\\|/g' ${z}) # change line breaks to '|'
	    grep ${linesmotu} ${input_file} > ${z}_grep&
        
done
#wait
echo "before cat"
cat ${temporal_dir}*_grep >> ${motus_tab_dir}${var}
echo "before MOTU"
echo "MOTU ${i} finished" 
rm -rf ${temporal_dir}