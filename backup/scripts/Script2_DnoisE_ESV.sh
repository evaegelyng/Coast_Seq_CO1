#!/bin/bash


# set many variables manually depending on your computer directories and files
vsearch_file=/home/evaes/eDNA/faststorage/Velux/CoastSequence/Spring/LerayXT/MJOLNIR/tmp/COSQ_vsearch.fasta # fasta file before swarm
tmp_dir=/home/evaes/eDNA/faststorage/Velux/CoastSequence/Spring/LerayXT/MJOLNIR/tmp/ #EES: directory where OTU table is stored
input_dir=/home/evaes/eDNA/faststorage/Velux/CoastSequence/Spring/LerayXT/MJOLNIR/results/ # directory where the output of swarm is stored (see line 46)
scripts_dir=/home/evaes/eDNA/faststorage/Velux/CoastSequence/Spring/LerayXT/MJOLNIR/scripts/  # directory where store required scripts (see line 163)
DnoisE_dir=/home/evaes/miniconda3/pkgs/dnoise-1.0-py38_0/lib/python3.8/site-packages/src/ 
cores=54 # number of cores (remember that there is a trade-off between speed and RAM. To much speed could kill the process due to RAM

MOTUS2RUN=${input_dir}COSQ_non_singleton_motu_list.txt # list of MOTUs



# create some directories to store all files separately
motus_dir=/home/evaes/eDNA/faststorage/Velux/CoastSequence/Spring/LerayXT/MJOLNIR/tmp/motus/
motus_tab_dir=/home/evaes/eDNA/faststorage/Velux/CoastSequence/Spring/LerayXT/MJOLNIR/tmp/motu_tab/
output_dir=/home/evaes/eDNA/faststorage/Velux/CoastSequence/Spring/LerayXT/MJOLNIR/tmp/output_Ad_corr/
final_dir=/home/evaes/eDNA/faststorage/Velux/CoastSequence/Spring/LerayXT/MJOLNIR/tmp/finalfiles/

if [[ ! -d ${motus_dir} ]]
 then
  mkdir ${motus_dir}
 fi

if [[ ! -d ${motus_tab_dir} ]]
 then
  mkdir ${motus_tab_dir}
 fi
 
if [[ ! -d ${output_dir} ]]
 then
  mkdir ${output_dir}
 fi

if [[ ! -d ${final_dir} ]]
 then
  mkdir ${final_dir}
 fi

RUN=false
if ${RUN}
 then
	# Size tag of the output of SWARM must be removed
	sed -i -e "s/;size=[0-9]*;/;/g" ${input_dir}COSQ_SWARM_output
	echo swarm_output modified

	# all names of MOTUs must be in the ${MOTUS2RUN} file

	# create a file per MOTU with all sequences that clustered into

	lines=$(wc -l ${MOTUS2RUN} | cut -f1 -d ' ')
	for i in $(seq 1 ${lines}) 
	do 
  		var=$(sed "${i}q;d" ${MOTUS2RUN})
  		grep ${var} ${input_dir}COSQ_SWARM_output >${motus_dir}${var}
  		sed -i -e "s/; /\\n/g" ${motus_dir}${var}
  		sed -i -e "s/;//g" ${motus_dir}${var}
 	done
	wait
	
	echo motu composition done

fi

RUN=true
if ${RUN}
 then
	# generate a .tab file for each MOTU with all sample information
	lines=$(wc -l ${MOTUS2RUN} | cut -f1 -d ' ')
	for i in $(seq 1 ${lines})
	do
  		var=$(sed "${i}q;d" ${MOTUS2RUN})
  		sed "1q;d" ${tmp_dir}COSQ_new.tab >${motus_tab_dir}${var}
  		mkdir ${motus_tab_dir}${var}_dir
  		temporal_dir=${motus_tab_dir}${var}_dir/
  		cd ${temporal_dir}
  		split -l200 -d ${motus_dir}${var}
  		for z in *
   		do
    			linesmotu=$(sed -e ':a;N;$!ba;s/\n/\\|/g' ${temporal_dir}${z}) # change line breaks to '|'
    			grep ${linesmotu} ${tmp_dir}COSQ_new.tab >${temporal_dir}${z}_grep&
   		done
  		wait
  		cat ${temporal_dir}*_grep >>${motus_tab_dir}${var}
  		echo MOTU ${i} acabat | rm -r ${temporal_dir} & 
 	done
	echo now wait untill all .tab files are made
	wait
	echo motu_tabs finished

fi

# here the step to retrieve entropy values from the whole dataset
RUN=false
if ${RUN}
then
	python3 ${DnoisE_dir}DnoisE.py --fasta_input ${vsearch_file} --csv_output ${input_dir} -n size -g
fi

RUN=false
if ${RUN}
then

	# run DnoisE

	lines=$(wc -l ${MOTUS2RUN} | cut -f1 -d ' ')
	for i in $(seq 1 ${lines})
 	do 
  		var=$(sed "${i}q;d" ${MOTUS2RUN})
  		if [ -f "${output_dir}${var}_Adcorr_denoised_ratio_d.csv" ]
   		then 
    			echo "${output_dir}${var}_Adcorr_denoised_ratio_d.csv exist" &
   		else
    			if [ $i == 1 ]
     			then
      				python3 ${DnoisE_dir}DnoisE.py --csv_input ${motus_tab_dir}${var} --csv_output ${output_dir}${var} -c ${cores} -s 4 -z 78 -a 4 -n count -e 0.4812,0.2407,1.0285 -y 2>${output_dir}${var}_Adcorr_progress &
     			else
      				python3 ${DnoisE_dir}DnoisE.py --csv_input ${motus_tab_dir}${var} --csv_output ${output_dir}${var} -c ${cores} -s 4 -z 78 -a 4 -n count -e 0.4812,0.2407,1.0285 -y 2>${output_dir}${var}_Adcorr_progress &
    			fi
  		fi
 	done
	wait

	echo DnoisE finished

fi

RUN=false
if ${RUN}
then

	lines=$(wc -l ${MOTUS2RUN} | cut -f1 -d ' ')

	for i in $(seq 1 ${lines})
 	do
  		motu=$(sed "${i}q;d" ${MOTUS2RUN})
  		if [ $i == 1 ]
   		then
    			head -n 1 ${output_dir}${motu}_Adcorr_denoised_ratio_d.csv > ${output_dir}ESV_Adcorr.csv
    			sed -i 's/id/motu,id/g' ${output_dir}ESV_Adcorr.csv
   		fi
  		awk 'NR>1' ${output_dir}${motu}_Adcorr_denoised_ratio_d.csv | awk -v var="$motu," '{print var $0;}' >>${output_dir}ESV_Adcorr.csv
 	done
 
fi

RUN=false
if ${RUN}
then

	# apply filters
	# 1. same filter of minimum relative abundance applied to MOTUs
	# 2. remove sequences with different length than the modal for each MOTU, also only lengths of 313+-(n*3) are allowed
	# 3. remove numts, this is:
		# a. seqs with codon stops
		# b. for metazoan and 313bp length sequences, sequences with changes in 5 well preserved aa.
		# (the genetic code used among the mitochondrial is the one with less codon stops per read and for metazoans also the one with less aa changes in 5 well preserved positions)
	# 4. remove ESV with less than 5 reads


	echo start filtering of ESV

	Rscript ${scripts_dir}Script2.1_DnoisE_ESV.R -i ${output_dir}ESV_Adcorr.csv -o ${final_dir} -a ${input_dir}metadata.RData -b ${input_dir}COSQ.Curated_LULU_euk.csv &


	wait
	echo finished

fi

