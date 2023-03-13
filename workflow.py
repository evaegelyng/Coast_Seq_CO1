from gwf import Workflow
import os, sys
import math
from glob import glob
import io
import pandas as pd
import numpy as np

project_name = "COSQ"

gwf = Workflow(defaults={"account": "edna"}) 

# Following workflow_part1, a complete metadata file was produced
# using the script make_metadata.R. Then, the dataset was cleaned 
# based on blank controls using clean_up_ASV_wise.R, and by removing
# sequences found in a single PCR replicate with no_sing_ASV_wise.R.
# The cleaned OTU table was then merged with the taxonomy from the
# output file of MJOLNIR (COSQ_final_dataset) using merge_otu_table_w_classified.
# The taxonomic identifications were then manually checked for sequences
# with hits of 97% similarity or more. Finally, OTUs that could be 
# confidently identified to marine species, contained at least 2 ASVs
# and were found in at least 10 clusters were selected using select_OTUs.R
# These 143 OTUs were then cleaned with DnoisE with the current workflow.
# After this workflow, Script2.1_DnoisE_ESV_230310.r was used to further filter
# OTUs and remove NUMTs

# Generate a .tab file for each MOTU with all sample information

motus_dir = "tmp/motus"
motus_tab_dir="tmp/motu_tab"
selected_motus="results/{}_pident_97_selected.txt".format(project_name)

with open(selected_motus, 'r') as fp:
    read = fp.readlines() 
    lines = len(read)

CORES=8

for i in range(0,len(read)):
    motu = read[i].strip()
    input_file = "tmp/{}_new.tab".format(project_name)
    output_file = "{}/{}".format(motus_tab_dir,motu)
    log_file = "{}/{}.log".format(motus_tab_dir,motu)

    gwf.target(
                name="tab_{}".format(motu),
                inputs=input_file,
                outputs=[output_file,log_file],
                cores=CORES,
                memory="8g",
                walltime="16:00:00",
            ) << """
                mkdir -p {motus_tab_dir}
                # parallelization with parsed outputs, option -k orders the outputs the same way as the inputs
                cat {motus_dir}/{motu} | parallel -j {CORES} -k --compress "python ./scripts/tab2.py {{}} {input_file}" > {motus_tab_dir}/{motu}
                # add header from the tabular file with motu, print output in a temporary file
                sed "1q;d" {input_file} | cat - {motus_tab_dir}/{motu} > {motus_tab_dir}/{motu}.tmp
                # substitute the file without header with the temporary file created above
                mv -f {motus_tab_dir}/{motu}.tmp {motus_tab_dir}/{motu}
                # print the log file for this target of the pipeline
                echo "hello" > {motus_tab_dir}/{motu}.log
            """.format(CORES=CORES, motu=motu, motus_tab_dir=motus_tab_dir, i=i, input_file=input_file, motus_dir=motus_dir,log_file=log_file)

# Run DnoisE. Remember to input entropy values from previous target to dnoise.sh
DnoisE_dir = "/home/evaes/miniconda3/pkgs/dnoise-1.1-py38_0/lib/python3.8/site-packages/src" 

output_dir = "tmp/output_Ad_corr"

cores = 2

for i in range(0,len(read)):
    motu = read[i].strip()
    
    input_files = ["{}/{}".format(motus_tab_dir,motu),"{}/{}.log".format(motus_tab_dir,motu)]
    
    output_file = "{}/{}_Adcorr_denoised_ratio_d.csv".format(output_dir,motu)
        
    gwf.target(
                name=f"dnoise_{motu}",
                inputs=input_files,
                outputs=output_file,
                cores=2,
                memory="16g",
                walltime="12:00:00",
            ) << """
                eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
                conda activate dnoise3
                mkdir -p {output_dir}
                scripts/dnoise.sh {motu} {output_dir} {DnoisE_dir} {motus_tab_dir} {cores}
            """.format(motu=motu, output_dir=output_dir, DnoisE_dir=DnoisE_dir, motus_tab_dir=motus_tab_dir, cores=cores)
