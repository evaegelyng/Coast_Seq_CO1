#!/bin/bash

input_file=$1 # SWARM output

MOTUS2RUN=$2 # list of MOTUs

# create output directory
motus_dir=$3

# Size tag of the output of SWARM must be removed
sed -i -e "s/;size=[0-9]*;/;/g" ${input_file}
echo swarm_output modified

# all names of MOTUs must be in the ${MOTUS2RUN} file

# create a file per MOTU with all sequences that clustered into

lines=$(wc -l ${MOTUS2RUN} | cut -f1 -d ' ')
for i in $(seq 1 ${lines}) 
do 
    var=$(sed "${i}q;d" ${MOTUS2RUN}) #Should change "var" to "motu" for consistency with following scripts
	grep ${var} ${input_file} >${motus_dir}/${var}
	sed -i -e "s/; /\\n/g" ${motus_dir}/${var}
	sed -i -e "s/;//g" ${motus_dir}/${var}
done
wait
	
echo motu composition done