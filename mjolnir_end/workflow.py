from gwf import Workflow
import os, sys
import math
from glob import glob
import io

project_name = "COSQ"

gwf = Workflow(defaults={"account": "edna"}) 

#Reformat tag files to ngsfilter files (obitools format)

libraries = [x for x in glob("/faststorage/project/eDNA/Velux/CoastSequence/Spring/LerayXT/backup/data/raw_data/*") if os.path.isdir(x)]

for library_root in libraries:
    library_id = os.path.basename(library_root)
    input_file = library_root + "/tags.txt"
    meta_file = "tmp/{}/metadata_{}.tsv".format(library_id,library_id) 
    output_file = "tmp/{}/ngsfilter_{}.tsv".format(library_id,library_id)
    
    gwf.target(
                name="reformat_{}_{}".format(project_name, library_id),
                inputs=input_file,
                outputs=[meta_file, output_file],
                cores=1,
                memory="1g",
                walltime="00:01:00",
            ) << """
                eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
                conda activate mjolnir
                mkdir -p tmp/{library_id}
                Rscript ./scripts/tags_to_ngsfilter.r {input_file} {library_id} {meta_file} {output_file}
            """.format(input_file=input_file, library_id=library_id, meta_file=meta_file, output_file=output_file)
          
#Using Mjolnir pipeline (https://github.com/uit-metabarcoding/MJOLNIR) to split fastq files for parallel processing, do paired-end alignment, demultiplexing, and read-length filtering.

for library_root in libraries:
    library_id = os.path.basename(library_root)
    input_files = []
    
    input_files.append(library_root + "/{}_1.fq".format(library_id))
    input_files.append(library_root + "/{}_2.fq".format(library_id))
    input_files.append("tmp/{}/ngsfilter_{}.tsv".format(library_id,library_id))
    
    output_files = []
    
    output_files.append("tmp/{}/{}_R1_part_01.fastq".format(library_id,library_id))
    output_files.append("tmp/{}/{}.filtered_length_part_01.fasta".format(library_id,library_id))
    output_files.append("tmp/{}/{}_sample_002.fasta".format(library_id,library_id)) #Use sample_002 as indication of completion. Sample_001 is a negative control in most libraries and often gives no fasta file.
        
    gwf.target(
                name="mjolnir_{}_{}".format(project_name, library_id),
                inputs=input_files,
                outputs=output_files,
                cores=16,    
                memory="8g",
                walltime="8:00:00",
            ) << """
                eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
                conda activate mjolnir
                cd tmp/{library_id}
                Rscript ./scripts/mjolnir.r {R1} {library_id}
                [[ -f "{library_id}_sample_002.fasta" ]] || touch "{library_id}_sample_002.fasta"
            """.format(R1=input_files[0], library_id=library_id) 
 
#De-replicate sequences and remove chimeras

for library_root in libraries:
    library_id = os.path.basename(library_root) 
    input_files = "tmp/{}/{}_sample_002.fasta".format(library_id,library_id)

    output_files = []

    output_files.append("tmp/{}/{}.vsearch.fasta".format(library_id,library_id))
    output_files.append("tmp/{}/{}.no_chimeras.fasta".format(library_id,library_id))
    output_files.append("tmp/{}/{}.new.fasta".format(library_id,library_id))
    output_files.append("tmp/{}/{}.new.tab".format(library_id,library_id))
    output_files.append("tmp/{}/{}.unique.fasta".format(library_id,library_id))

    gwf.target(
                name="hela_{}_{}".format(project_name, library_id), 
                inputs=input_files,
                outputs=output_files,
                cores=3,
                memory="8g",
                walltime="2:00:00",
            ) << """
                eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
                conda activate mjolnir
                cd tmp/{library_id}
                Rscript ./scripts/hela.r {library_id}
            """.format(library_id=library_id)      

#Combine fasta files of unique non-chimera sequences from all libraries

input_files = []
for library_root in libraries:
    library_id = os.path.basename(library_root)  
    input_files.append("tmp/{}/{}.unique.fasta".format(library_id,library_id))

output_files = []
 
output_files.append("tmp/{}.no_chimeras.fasta".format(project_name))
output_files.append("tmp/{}.unique.fasta".format(project_name))
output_files.append("tmp/{}.new.fasta".format(project_name))
output_files.append("tmp/{}.vsearch.fasta".format(project_name))
output_files.append("tmp/{}.new.tab".format(project_name))

all_files = glob("/faststorage/project/eDNA/Velux/CoastSequence/Spring/LerayXT/MJOLNIR/tmp/*/????.unique.fasta")

old_stdout = sys.stdout
new_stdout = io.StringIO()
sys.stdout = new_stdout
print(" ".join(all_files)) 
output = new_stdout.getvalue()
sys.stdout = old_stdout   
output = output.rstrip("\n") 
    
gwf.target(
            name="hela_all_{}".format(project_name),
            inputs=input_files,
            outputs=output_files,
            cores=2,
            memory="196g",
            walltime="12:00:00",
        ) << """
            eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
            conda activate mjolnir
            
            cat {output} > tmp/{project_name}.no_chimeras.fasta
            obiuniq -m sample tmp/{project_name}.no_chimeras.fasta | obigrep -p 'count>1' > tmp/{project_name}.unique.fasta
            echo "HELA will change sequence identifiers to a short index"
            obiannotate --seq-rank tmp/{project_name}.unique.fasta | obiannotate --set-identifier \'\"\'"{project_name}"\'_%09d\" % seq_rank\' > tmp/{project_name}.new.fasta
            echo "HELA will change the format to vsearch, so ODIN can use it for SWARM."
            Rscript ./scripts/obi2vsearch.r tmp/{project_name}
            echo "File tmp/{project_name}.vsearch.fasta written."
            echo "HELA is obtaining a table file with abundances of unique sequence in each sample"
            obitab -o tmp/{project_name}.new.fasta >  tmp/{project_name}.new.tab
            echo "HELA is done."
        """.format(output=output, project_name=project_name) 

#Using Mjolnir pipeline to perform OTU clustering and denoising. Run separately, as DnoisE was not compatible with the "mjolnir" conda environment. The SWARM clustering itself was done by Owen Wangensteen at UiT, producing the four commented output files below. Replaced "." with "_" in input file names.
 
input_files = []

input_files.append("tmp/{}_new.tab".format(project_name))
input_files.append("tmp/{}_vsearch.fasta".format(project_name))

output_files = []
 
#output_files.append("results/{}_SWARM_seeds.fasta".format(project_name))
#output_files.append("results/{}_SWARM13nc_stats".format(project_name))
#output_files.append("results/{}_SWARM_output".format(project_name))
#output_files.append("results/{}_non_singleton_motu_list.txt".format(project_name))
output_files.append("results/{}_SWARM_output.counts.csv".format(project_name))
output_files.append("results/{}_SWARM_output.ESV.csv".format(project_name))

gwf.target(
            name="odin_{}".format(project_name),
            inputs=input_files,
            outputs=output_files,
            cores=18,
            memory="384g",
            walltime="7-00:00:00",
        ) << """
            eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
            conda activate dnoise_2
            cd results
            Rscript ../scripts/odin_220224.r ../tmp/{project_name} {project_name}
        """.format(project_name=project_name) 
            
#Using obigrep to remove singletons (from last part of ODIN function). 

input_file = "results/{}_SWARM_seeds.fasta".format(project_name)
output_file = "results/{}_seeds_abundant.fasta".format(project_name)

gwf.target(
            name="nonsingleton_{}".format(project_name),
            inputs=input_file,
            outputs=output_file,
            cores=1,
            memory="1g",
            walltime="00:01:00",
        ) << """
            eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
            conda activate mjolnir
            echo "ODIN will now remove MOTUs with total abundance less than ",min_reads_MOTU," from the fasta output file, to decrease THOR's workload."
            sed -i 's/;size/ size/g' {input_file}
            obigrep -p 'size>1' {input_file} > {output_file}
            echo "ODIN is done."
        """.format(input_file=input_file, output_file=output_file) 
        
#Assign the taxonomy to the representative sequence of each MOTU. 

input_file = "results/{}_seeds_abundant.fasta".format(project_name)
output_file = "results/{}_ecotag_annotated.tsv".format(project_name)

gwf.target(
            name="thor_{}".format(project_name),
            inputs=input_file,
            outputs=output_file,
            cores=54,
            memory="196g",
            walltime="3-00:00:00",
        ) << """
            eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
            conda activate mjolnir
            cd results
            Rscript ../scripts/thor.r {project_name}
        """.format(project_name=project_name)         
