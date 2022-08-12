#!/bin/bash

motus_tab_dir=$1
i=$2
MOTUS2RUN=$3
input_file=$4
motus_dir=$5

var=$(sed "${i}q;d" ${MOTUS2RUN}) #get name of MOTU, corresponding to the name of representative sequence 
sed "1q;d" ${input_file} > ${motus_tab_dir}/${var} #get all sample names from the tab file of all MOTUs
mkdir -p ${motus_tab_dir}/${var}_dir #make a temporary folder for this MOTU
temporal_dir=${motus_tab_dir}/${var}_dir/
cd ${temporal_dir} #go into the temporary folder and split the file of sequences (sequence names) belonging to this MOTU
split -l200 -d ../../../${motus_dir}/${var}
cd - #go back to root folder
for z in `ls ${temporal_dir}/*`
do
        
	    linesmotu=$(sed -e ':a;N;$!ba;s/\n/\\|/g' ${z}) # change line breaks to '|'
	    echo $linesmotu | xargs -n 1 -d "|" | xargs -0 grep {} $input_file > ${z}_grep & #extract read counts for this sequence name (z) and run in background
             
done
wait #wait for all the grep commands to be done
echo "before cat"
cat ${temporal_dir}*_grep >> ${motus_tab_dir}/${var}
echo "before MOTU"
echo "MOTU ${i} finished" 
rm -rf ${temporal_dir}