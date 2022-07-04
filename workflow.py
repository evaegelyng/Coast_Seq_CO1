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
                name="ran_{}_{}".format(project_name, library_id),
                inputs=input_files,
                outputs=output_files,
                cores=16,    
                memory="8g",
                walltime="8:00:00",
            ) << """
                eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
                conda activate mjolnir
                cd tmp/{library_id}
                Rscript ./scripts/ran_freyja.r {R1} {library_id}
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
output_files.append("tmp/{}_vsearch.fasta".format(project_name))
output_files.append("tmp/{}_new.tab".format(project_name))

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
 
# SWARM was run separately by O. Wangensteen at UiT, producing the following outputs. However, it could also be run on GenomeDK - this just require 72 cores and almost 2 weeks..
#output_files.append("results/{}_SWARM_seeds.fasta".format(project_name))
#output_files.append("results/{}_SWARM13nc_stats".format(project_name))
#output_files.append("results/{}_SWARM_output".format(project_name))
#output_files.append("results/{}_non_singleton_motu_list.txt".format(project_name))

# DnoisE is run in a separate target for now, producing the following output
#output_files.append("results/{}_SWARM_output.ESV.csv".format(project_name))

output_files.append("results/{}_SWARM_output_counts.tsv".format(project_name))

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

# Use DnoisE to remove likely erroneous sequences

input_files = []

input_files.append("tmp/{}_vsearch.fasta".format(project_name))
input_files.append("results/{}_SWARM_output".format(project_name))

output_files = []

output_files.append("results/{}_SWARM_output.ESV.csv".format(project_name))

gwf.target(
            name="dnoise_{}".format(project_name),
            inputs=input_files,
            outputs=output_files,
            cores=54,
            memory="196g",
            walltime="12:00:00",            
        ) << """
            eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
            conda activate dnoise_2
            cd results
            ../scripts/Script2_DnoisE_ESV.sh
        """.format(project_name=project_name) 
            
#Using obigrep to remove singletons (from last part of ODIN function). 

input_file = "results/{}_SWARM_seeds.fasta".format(project_name)
output_file = "results/{}_seeds_abundant.fasta".format(project_name)

gwf.target(
            name="nonsingleton_{}_{}".format(project_name, library_id),
            inputs=input_file,
            outputs=output_file,
            cores=8,
            memory="8g",
            walltime="2:00:00",
        ) << """
            eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
            conda activate mjolnir
            echo "ODIN will now remove MOTUs with total abundance less than ",min_reads_MOTU," from the fasta output file, to decrease THOR's workload."
            obigrep -p 'size>1' {input_file} > {output_file}
            echo "ODIN is done."
        """.format(library_id=library_id, input_file=input_file, output_file=output_file) 

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

#We decided to also perform a BLAST search against our own database of eukaryote CO1 sequences from BOLD and GenBank, and taxonomic assignment using our own R script (see MetaBarFlow pipeline on GitHub).

###Split fasta file (the nochim one with chimeras removed) into K parts
def splitter(inputFile, K=99):
    inputs = [inputFile]
    outputs = ["tmp/split/split.log.txt"]
    options = {
        'cores': 1,
        'memory': '2g',
        'walltime': '1:00:00'
    }
    spec = '''
    eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
    conda activate metabar_2021
    seqkit split -O tmp/split/ {inputFile} -p {K} -2
    echo "hello" > tmp/split/split.log.txt
    '''.format(inputFile=inputFile, K=K)
    return inputs, outputs, options, spec

#####blast a single k-th file
def blaster(k, outFolder):
    inputFasta = 'tmp/split/{}_seeds_abundant.part_'.format(project_name)+'{:0>3d}'.format(k)+'.fasta'
    inputs = [inputFasta]
    outBlast = outFolder + '/blast.' + str(k) + '.blasthits'
    outLog = outFolder + '/blast.' + str(k) + '.txt'
    outputs = [
      outBlast,
      outLog
    ]
    options = {
        'cores': 2,
        'memory': '32g',
        'walltime': '12:00:00'
    }
    spec = '''
    eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
    conda activate metabar_2021
    export BLASTDB=/faststorage/project/eDNA/blastdb/Eukaryota_COI_NOBAR/BLAST_db/
    mkdir -p {out}
    echo "RUNNING THREAD {k} BLAST"
    blastn -db /faststorage/project/eDNA/blastdb/Eukaryota_COI_NOBAR/BLAST_db/Eukaryota_NOBAR.db -max_target_seqs 500 -num_threads 4 -outfmt "6 std qlen qcovs staxid" -out {outBlast} -qcov_hsp_perc 90 -perc_identity 80 -query {inputFasta}
    echo "hello" > {outLog}
    echo "DONE THREAD {k}"
    '''.format(out=outFolder, k=k, inputFasta=inputFasta, outBlast=outBlast, outLog=outLog)
    return inputs, outputs, options, spec

def taxonomy(taxonomyFolder, blastFolder, k):
    inputFile = blastFolder + '/blast.' + str(k) + '.blasthits'
    inputs = [inputFile , blastFolder + '/blast.' + str(k) + '.txt']
    summaryFile = taxonomyFolder + '/summary.' + str(k) + '.txt'
    outputFile = taxonomyFolder + '/taxonomy.' + str(k) + '.txt'
    outputs = [summaryFile, outputFile]
    options = {
        'cores': 1,
        'memory': '32g',
        'walltime': '4:00:00'
    }
    
    spec = '''
    eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
    conda activate metabar_2021
    mkdir -p {taxonomyFolder}
    # Check if blast file is empty
    if [ `cat {inputFile} | wc -l` != 0 ]
    then
      Rscript scripts/taxonomy_bold_nt_220601.r {inputFile} {summaryFile} {outputFile}
    else
      touch {outputFile}
    fi
    if grep -q "Query coverage is less than 100% for all hits" ".gwf/logs/taxonomy_" + str(k) + ".stdout"
    then
      touch {outputFile}
    if grep -q "Sequence identity is less than 90% for all hits" ".gwf/logs/taxonomy_" + str(k) + ".stdout"
    then
      touch {outputFile}
    '''.format(taxonomyFolder=taxonomyFolder, inputFile=inputFile, summaryFile=summaryFile, outputFile=outputFile) 
    return inputs, outputs, options, spec

inputName = "results/{}_seeds_abundant.fasta".format(project_name)

gwf.target_from_template( 'split', splitter(inputFile=inputName) )

parts=glob('tmp/split/{}_seeds_abundant.part*.fasta'.format(project_name))
K=len(parts)
                                                                
for k in range(1,K+1):
  gwf.target_from_template( 'blaster_{}'.format(k), blaster(k=k, outFolder='tmp/blast') )
  gwf.target_from_template( 'taxonomy_{}'.format(k), taxonomy(taxonomyFolder='tmp/taxonomy', blastFolder='tmp/blast', k=k) )

### Combine all the small taxonomical summary files into one large file

input_files = glob('tmp/taxonomy/summary*.txt')
 
output_files = ['results/summary.txt']
    
gwf.target(
   name="combine_summary_{}".format(project_name),
   inputs=input_files,
   outputs=output_files,
   cores=1,
   memory="1g",
   walltime="00:10:00",
 ) << """
    head -n1 tmp/taxonomy/summary.1.txt > results/summary.txt
    for fname in tmp/taxonomy/summary*.txt
    do
        tail -n +2 $fname >> results/summary.txt
    done
   """    
      
### Combine all the small taxonomical classfication files into one large file

input_files = glob('tmp/taxonomy/taxonomy*.txt')
 
output_files = ['results/classified.txt']
    
gwf.target(
   name="combine_taxonomy_{}".format(project_name),
   inputs=input_files,
   outputs=output_files,
   cores=1,
   memory="1g",
   walltime="00:10:00",
 ) << """
    head -n1 tmp/taxonomy/taxonomy.1.txt > results/classified.txt
    for fname in tmp/taxonomy/taxonomy*.txt
    do
        tail -n +2 $fname >> results/classified.txt
    done
   """         
   
# FRIGGA will integrate the information of MOTU abundances and taxonomy assignment from ODIN & THOR in a single table

input_files = []

input_files.append("results/{}_classified.tsv".format(project_name))
input_files.append("results/{}_SWARM_output_counts.tsv".format(project_name))
output_file = "results/{}_All_MOTUs_classified.tsv".format(project_name)

gwf.target(
            name="frigga_{}".format(project_name),
            inputs=input_files,
            outputs=output_file,
            cores=1,
            memory="8g",
            walltime="00:10:00",
        ) << """
            eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
            conda activate mjolnir
            cd results
            Rscript ../scripts/frigga.r {project_name}
        """.format(project_name=project_name)         
        
# LOKI will remove the pseudogenes and will keep track of the taxonomic information of the removed MOTUs
        
input_file= "results/{}_All_MOTUs_classified.tsv".format(project_name)
output_files = []

output_files.append("results/{}_match_list.txt".format(project_name))
output_files.append("results/{}_Discarded_LULU.tsv".format(project_name))
output_files.append("results/{}_Curated_LULU.tsv".format(project_name))
output_files.append("results/{}_Deleted_LULU_fate.tsv".format(project_name))

gwf.target(
            name="loki_{}".format(project_name),
            inputs=input_file,
            outputs=output_files,
            cores=1,
            memory="56g",
            walltime="12:00:00",
        ) << """
            eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
            conda activate mjolnir
            cd results
            Rscript ../scripts/loki.r {project_name}
        """.format(project_name=project_name)         
                        
### Combine all the metadata files into one large file

#First, copy all the input files to the tmp folder to allow using the glob function below. Could also just output directly to tmp when creating the files, but it seems more organized to have them in their corresponding library folder
#cd tmp
#for FOLDER_NAME in ./S*
#do
#  cp "$FOLDER_NAME"/metadata* . 
#  echo  "$FOLDER_NAME done"
#done

#for FOLDER_NAME in ./W*
#do
#  cp "$FOLDER_NAME"/metadata* . 
#  echo  "$FOLDER_NAME done"
#done

input_files = glob('tmp/metadata*')
 
output_file = "results/{}_metadata.tsv".format(project_name)
    
gwf.target(
   name="metadata_{}".format(project_name),
   inputs=input_files,
   outputs=output_file,
   cores=1,
   memory="8g",
   walltime="00:10:00",
 ) << """
    head -n1 tmp/metadata_S111.tsv > results/COSQ_metadata.tsv
    for fname in tmp/metadata*
    do
        tail -n +2 $fname >> results/COSQ_metadata.tsv
    done
   """   

#Now the copies of the metadata files in the tmp folder can be deleted
#cd tmp
#rm_force metadata*

# RAGNAROC will change the names of the samples to recover the original names and will remove unnecessary columns

input_files= []

input_files.append("results/{}_metadata.tsv".format(project_name))
input_files.append("results/{}_Curated_LULU.tsv".format(project_name))
input_files.append("results/{}_All_MOTUs_classified.tsv".format(project_name))

output_file = "results/{}_final_dataset_classified.tsv".format(project_name)

gwf.target(
            name="ragnaroc_{}".format(project_name),
            inputs=input_files,
            outputs=output_file,
            cores=1,
            memory="56g",
            walltime="4:00:00",
        ) << """
            eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
            conda activate mjolnir
            cd results
            Rscript ../scripts/ragnaroc.r {project_name}
        """.format(project_name=project_name)        
