#!/bin/bash

motu=$1
input_file=$2
output_file=$3
cores=$4

#motu=$(sed "${i}q;d" ${MOTUS2RUN}) #get name of MOTU, corresponding to the name of representative sequence 
sed "1q;d" ${input_file} > ${output_file} #get all sample names from the tab file of all MOTUs
#mkdir -p ${motus_tab_dir}/${motu}_dir #make a temporary folder for this MOTU
#temporal_dir=${motus_tab_dir}/${motu}_dir/
#cd ${temporal_dir} #go into the temporary folder and split the file of sequences (sequence names) belonging to this MOTU
#split -l200 -d ../../../${motus_dir}/${motu}
#cd - #go back to root folder
#rm -f ${temporal_dir}/*_grep

#for z in `ls ${temporal_dir}/*`
#do
	    #linesmotu=$(sed -e ':a;N;$!ba;s/\n/\\|/g' ${z}) # change line breaks to '|'
	    #echo $linesmotu | xargs -n 1 -d "|" | xargs -I {} grep {} $input_file > ${z}_grep & #extract read counts for this sequence name (z) and run in background
        #cat $z | xargs -P${cores} -I {} fgrep -m 1 {} $input_file > ${z}_grep
		#echo "Running with sequence $z" 
		#time python ./scripts/tab2.py ${z} ${temporal_dir} ${input_file}  
#done

cat ${motus_dir}/${motu} | xargs -P${cores} -I {} ./scripts/tab2.py {} ${input_file} >> ${motus_tab_dir}/${motu}

#wait #wait for all the grep commands to be done
#echo "merging all PCR counts"
#cat ${temporal_dir}/*_grep >> ${motus_tab_dir}/${motu}
echo "MOTU ${i} finished" 
#rm -rf ${temporal_dir}