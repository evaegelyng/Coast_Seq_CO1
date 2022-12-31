from gwf import Workflow
import os, sys
import math
from glob import glob
import io
import pandas as pd
import numpy as np

project_name = "COSQ"

gwf = Workflow(defaults={"account": "edna"}) 

#Using Mjolnir pipeline to perform OTU clustering.

db.total="results/db.total.tsv"

#read MOTUS from motus file
MOTUS = np.ravel(np.array(pd.read_csv(db.total, sep='\t')["id"]))
#batches of MOTU names in a dictionary
Lmotus = len(MOTUS) #nr of MOTUs

input_files = []

input_files.append("../tmp/{}_new.tab".format(project_name))
input_files.append("../tmp/{}_SWARM_output".format(project_name))

output_file = "results/{}_SWARM_output_counts.tsv".format(project_name)

gwf.target(
            name="odin_{}".format(project_name),
            inputs=input_files,
            outputs=output_file,
            cores=32,
            memory="364g",
            walltime="7-00:00:00",            
        ) << """
            eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
            conda activate dnoise_2
            cd results
            Rscript ../scripts/odin_test.r ../tmp/{project_name} {project_name}
        """.format(project_name=project_name) 
