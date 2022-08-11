#!/bin/bash

motus_tab_dir=$1
i=$2
MOTUS2RUN=$3
input_file=$4
motus_dir=$5

var=$(sed "${i}q;d" ${MOTUS2RUN}) #name COSQxxxx where xxxx is some number
sed "1q;d" ${input_file} > ${motus_tab_dir}/${var} #get the COSQxxxx names from the input file
mkdir -p ${motus_tab_dir}/${var}_dir #make a COSQxxxx_dir folder for each COSQ name
temporal_dir=${motus_tab_dir}/${var}_dir/
cd ${temporal_dir} #go into the temporal dir and split the COSQ file
split -l200 -d ../../../${motus_dir}/${var}
cd - #go back to previous folder
for z in `ls ${temporal_dir}/*`
do
        
	    linesmotu=$(sed -e ':a;N;$!ba;s/\n/\\|/g' ${z}) # change line breaks to '|'
	    grep ${linesmotu} ${input_file} > ${z}_grep
        
done
#wait
echo "before cat"
cat ${temporal_dir}*_grep >> ${motus_tab_dir}/${var}
echo "before MOTU"
echo "MOTU ${i} finished" 
rm -rf ${temporal_dir}