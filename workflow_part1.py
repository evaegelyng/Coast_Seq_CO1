from gwf import Workflow
import os, sys
import math
from glob import glob
import io
import pandas as pd
import numpy as np

project_name = "COSQ"

gwf = Workflow(defaults={"account": "edna"}) 

#Reformat tag files to ngsfilter files (obitools format)

libraries = [x for x in glob("/faststorage/project/eDNA/Velux/CoastSequence/Spring/LerayXT/backup/data/raw_data/*") if os.path.isdir(x)]

for library_root in libraries:
    library_id = os.path.basename(library_root)
    input_file = library_root + "/tags.txt"
    meta_file = "tmp/metadata_{}.tsv".format(library_id) 
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
#To determine a minimum read count threshold, the counts are extracted from the output of obiuniq

input_files = []
for library_root in libraries:
    library_id = os.path.basename(library_root)  
    input_files.append("tmp/{}/{}.unique.fasta".format(library_id,library_id))

output_files = []
 
output_files.append("tmp/{}.no_chimeras.fasta".format(project_name))
output_files.append("tmp/{}.unique.fasta".format(project_name))
output_files.append("tmp/{}.new.fasta".format(project_name))
output_files.append("tmp/{}.vsearch.fasta".format(project_name))
output_files.append("tmp/{}_new.tab".format(project_name))
output_files.append("tmp/{}_new.tab.names".format(project_name))

cores=2

gwf.target(
            name="hela_all_{}".format(project_name),
            inputs=input_files,
            outputs=output_files,
            cores=cores,
            memory="32g",
            walltime="12:00:00",
        ) << """
            eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
            conda activate mjolnir
            
            #concatenate into main fasta file
            echo "Concatenating " `find /faststorage/project/eDNA/Velux/CoastSequence/Spring/LerayXT/MJOLNIR/tmp/ -name "????.unique.fasta"|wc -l` " Files"
            #cat `find /faststorage/project/eDNA/Velux/CoastSequence/Spring/LerayXT/MJOLNIR/tmp/ -name "????.unique.fasta"` > tmp/{project_name}.no_chimeras.fasta
                        
            echo "Making unique values from the following fasta file " `ls -sh tmp/{project_name}.no_chimeras.fasta`
            #obiuniq -m sample tmp/{project_name}.no_chimeras.fasta > tmp/{project_name}.no_chimeras.unique.fasta
            
            #writing occurrences of each number of sequences counts
            #grep "^>" tmp/{project_name}.no_chimeras.unique.fasta | cut -d ";" -f2 > tmp/{project_name}.no_chimeras.unique.counts
            #sed -i 's/count=//g' tmp/{project_name}.no_chimeras.unique.counts
            #sort -n -S 50% --parallel={cores}  tmp/{project_name}.no_chimeras.unique.counts | uniq -c > tmp/{project_name}.no_chimeras.unique.occurrences.txt
            
            #finding the threshold using occurrences of each count value
            Rvalues=(`Rscript --vanilla scripts/threshold_choice.R tmp/{project_name}.no_chimeras.unique.occurrences.txt`)
            threshold=${{Rvalues[0]}}
            keep=${{Rvalues[1]}}
            totalcounts=${{Rvalues[2]}}

            echo "Keeping" $keep "sequences with more than" $threshold "counts out of" $totalcounts "sequences"
            parameter=`echo "'count>"$threshold"'"`
            cmd="obigrep -p $parameter tmp/{project_name}.no_chimeras.unique.fasta > tmp/{project_name}.unique.fasta"
            #eval $cmd
            
            echo "HELA will change sequence identifiers to a short index"
            obiannotate --seq-rank tmp/{project_name}.unique.fasta | obiannotate --set-identifier \'\"\'"{project_name}"\'_%09d\" % seq_rank\' > tmp/{project_name}.new.fasta
            
            echo "HELA will change the format to vsearch, so ODIN can use it for SWARM."
            Rscript ./scripts/obi2vsearch.r tmp/{project_name}
            
            echo "File tmp/{project_name}.vsearch.fasta written."
            echo "HELA is obtaining a table file with abundances of unique sequence in each sample"
            obitab -o tmp/{project_name}.new.fasta >  tmp/{project_name}_new.tab
            
            echo "HELA is done."
            echo "Making a file containing only sequence names, used for indexing in tab2.py script"
            cut -f1 tmp/{project_name}_new.tab > tmp/{project_name}_new.tab.names
        """.format(project_name=project_name, cores=cores) 

#Using Mjolnir pipeline to perform OTU clustering.
#NB! Move last line of R script to workflow
 
input_files = []

input_files.append("tmp/{}_new.tab".format(project_name))
input_files.append("tmp/{}.vsearch.fasta".format(project_name))

output_files = []
 
output_files.append("results/{}_SWARM_seeds.fasta".format(project_name))
output_files.append("results/{}_SWARM13nc_stats".format(project_name))
output_files.append("results/{}_SWARM_output".format(project_name))
output_files.append("results/{}_non_singleton_motu_list.txt".format(project_name))
output_files.append("results/{}_SWARM_output_counts.tsv".format(project_name))

gwf.target(
            name="odin_{}".format(project_name),
            inputs=input_files,
            outputs=output_files,
            cores=32,
            memory="364g",
            walltime="7-00:00:00",            
        ) << """
            eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
            conda activate dnoise_2
            cd results
            Rscript ../scripts/odin_220224.r ../tmp/{project_name} {project_name}
        """.format(project_name=project_name) 

# Count number of reads per MOTU
#grep "^>" results/COSQ_SWARM_seeds.fasta | cut -d " " -f2 > results/COSQ_SWARM_seeds.counts
#sed -i 's/size=//g' results/COSQ_SWARM_seeds.counts
#sed -i 's/\;//g' results/COSQ_SWARM_seeds.counts
#sort -n -S 50% --parallel=1  results/COSQ_SWARM_seeds.counts | uniq -c > results/COSQ_SWARM_seeds.counts.occurences.txt

# Create a file per MOTU with all sequences that clustered into. 

input_file = "results/{}_SWARM_output".format(project_name)
cleanfile="results/{}_final_dataset_cleaned_pident_70.tsv".format(project_name) 
motus2run_file="results/{}_MOTUS2RUN.csv".format(project_name)
MOTUS2RUN = pd.read_csv(cleanfile, sep='\t', quoting=False)["id"]
MOTUS2RUN.to_csv(motus2run_file, index=False, header=False)

motus_dir = "tmp/motus"
output_file = "{}/{}_000000015".format(motus_dir,project_name) # Only added the first MOTU file as an output. Check that log file says "motu composition done"

gwf.target(
            name="motus_{}".format(project_name),
            inputs=[input_file,cleanfile],
            outputs=output_file,
            cores=1,
            memory="2g",
            walltime="04:00:00",            
        ) << """
            mkdir -p {motus_dir}
            scripts/motus.sh {input_file} {motus2run_file} {motus_dir}
        """.format(project_name=project_name, input_file=input_file, output_file=output_file, motus2run_file=motus2run_file, motus_dir=motus_dir) 

#selected_motus="results/{}_final_dataset_cleaned_pident_97_selected.tsv".format(project_name) #Need to replace with final file
#read MOTUS from motus file
#MOTUS = np.ravel(np.array(pd.read_csv(selected_motus, sep='\t', quoting=False)["id"]))
#shuffle MOTUS
#np.random.seed(123)
#np.random.shuffle(MOTUS)
#batches of MOTU names in a dictionary
#Lmotus = len(MOTUS) #nr of MOTUs
#N = 2 #nr of MOTUs in a batch
#B = Lmotus // N #nr batches
#B_last = Lmotus % N #size of last batch
#MOTUSfiles = [] #empty string list with B or B+1 batches

#add a string of N motu paths in each entry
#for i in range(B): 
#    vec = [ f"{MOTUS[i*N+j]}" for j in range(N) ]
#    vec2 = ' '.join(vec)
#    MOTUSfiles.append(vec2)
#add the last B_last motu if necessary
#if B_last>0:
#    vec = [ f"{MOTUS[B*N+j]}" for j in range(B_last) ]
#    vec2 = ' '.join(vec)
#    MOTUSfiles.append(vec2)

#CORES=32
#for batch in range(len(MOTUSfiles)):
#    MOTUSlist = MOTUSfiles[batch]
#    MOTUSlist = MOTUSlist.split(' ')


#    input_file = "tmp/{}_new.tab".format(project_name)
#    output_file = ["{}/BATCH_{}.log".format(motus_tab_dir,batch)] + ["{}/{}".format(motus_tab_dir,i) for i in MOTUSlist]
    

#    gwf.target( name=f"tab_batch_{batch}",
#                inputs=input_file,
#                outputs=output_file,
#                cores=CORES,
#                memory="8g",
#                walltime="7-00:00:00",
#            ) << """ 
#                mkdir -p {motus_tab_dir}
#                rm -f output_file
#                for NAME in {MOTUSlist}
#                do
#                    FILENAME=$(echo $NAME | sed 's/,//g' | sed 's/[][]//g')
#                    if [ ! -f {motus_tab_dir}/${{FILENAME}}.log ]
#                    then    
                        # parallelization with parsed outputs, option -k orders the outputs the same way as the inputs
#                        cat {motus_dir}/$FILENAME | parallel -j {CORES} -k --compress "python ./scripts/tab2.py {{}} {input_file}" > {motus_tab_dir}/$FILENAME
                        # add header from the tabular file with motu, print output in a temporary file
#                        sed "1q;d" {input_file} | cat - {motus_tab_dir}/$FILENAME > {motus_tab_dir}/$FILENAME.tmp
                        # substitute the file without header with the temporary file created above
#                        mv -f {motus_tab_dir}/$FILENAME.tmp {motus_tab_dir}/$FILENAME
                        # print the log file for this target of the pipeline
#                        echo "hello" > {motus_tab_dir}/${{FILENAME}}.log
#                    fi
#                done
#                echo "hello" > {motus_tab_dir}/BATCH_{batch}.log
#            """.format(CORES=CORES, MOTUSlist=MOTUSlist,
#            motus_tab_dir=motus_tab_dir, motus_dir=motus_dir,
#            input_file=input_file, batch=batch)

# Remove trailing spaces from COSQ_vsearch.fasta

input_file = "tmp/{}.vsearch.fasta".format(project_name)

output_file = "tmp/{}.vsearch_no_spaces.fasta".format(project_name)

gwf.target(
            name="spaces_{}".format(project_name),
            inputs=input_file,
            outputs=output_file,
            cores=1,
            memory="4g",
            walltime="00:10:00",            
        ) << """
            scripts/remove_spaces.sh {input_file} {output_file}
        """.format(project_name=project_name, input_file=input_file, output_file=output_file) 

# here the step to retrieve entropy values from the whole dataset

DnoisE_dir = "/home/evaes/miniconda3/pkgs/dnoise-1.1-py38_0/lib/python3.8/site-packages/src" 

vsearch_file = "tmp/{}.vsearch_no_spaces.fasta".format(project_name)

output_dir = "results"

output_file = "results/results_entropy_values.csv"

gwf.target(
            name="entropy_{}".format(project_name),
            inputs=vsearch_file,
            outputs=output_file,
            cores=1,
            memory="196g",
            walltime="12:00:00",            
        ) << """
            eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
            conda activate dnoise3
            python3 {DnoisE_dir}/DnoisE.py --fasta_input {vsearch_file} --csv_output {output_dir} -n size -g
            mv results_entropy_values.csv results
        """.format(project_name=project_name, DnoisE_dir=DnoisE_dir, vsearch_file=vsearch_file, output_dir=output_dir)

###Split fasta file into K parts
#NB! Should add removal of old index file, as this seems to not be overwritten
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
    inputFasta = 'tmp/split/{}_SWARM_seeds.part_'.format(project_name)+'{:0>3d}'.format(k)+'.fasta'
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
    export BLASTDB=/faststorage/project/eDNA/blastdb/Euk_COI_NOBAR_2022/BLAST_db/
    mkdir -p {out}
    echo "RUNNING THREAD {k} BLAST"
    blastn -db /faststorage/project/eDNA/blastdb/Euk_COI_NOBAR_2022/BLAST_db/Eukaryota_NOBAR.db -max_target_seqs 500 -num_threads 4 -outfmt "6 std qlen qcovs staxid" -out {outBlast} -qcov_hsp_perc 90 -perc_identity 70 -query {inputFasta}
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
      Rscript scripts/taxonomy_v0.1.1_pident_80.r {inputFile} {summaryFile} {outputFile}
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

inputName = "results/{}_SWARM_seeds.fasta".format(project_name)

gwf.target_from_template( 'split', splitter(inputFile=inputName) )

parts=glob('tmp/split/{}_SWARM_seeds.part*.fasta'.format(project_name))
K=len(parts)
                                                                
for k in range(1,K+1):
  gwf.target_from_template( 'blaster_{}'.format(k), blaster(k=k, outFolder='tmp/blast') )
  gwf.target_from_template( 'taxonomy_{}'.format(k), taxonomy(taxonomyFolder='tmp/taxonomy', blastFolder='tmp/blast', k=k) )

### Combine all the small taxonomical summary files into one large file

input_files = glob('tmp/taxonomy/summary*.txt')
 
output_file = 'results/{}_summary.txt'.format(project_name)
    
gwf.target(
   name="combine_summary_{}".format(project_name),
   inputs=input_files,
   outputs=output_file,
   cores=1,
   memory="1g",
   walltime="00:10:00",
 ) << """
    head -n1 tmp/taxonomy/summary.1.txt > results/{project_name}_summary.txt
    for fname in tmp/taxonomy/summary*.txt
    do
        tail -n +2 $fname >> results/{project_name}_summary.txt
    done
   """.format(project_name=project_name)    
      
### Combine all the small taxonomical classfication files into one large file

input_files = glob('tmp/taxonomy/taxonomy*.txt')
 
output_file = 'results/{}_classified.tsv'.format(project_name)
    
gwf.target(
   name="combine_taxonomy_{}".format(project_name),
   inputs=input_files,
   outputs=output_file,
   cores=1,
   memory="1g",
   walltime="00:10:00",
 ) << """
    head -n1 tmp/taxonomy/taxonomy.1.txt > results/{project_name}_classified.tsv
    for fname in tmp/taxonomy/taxonomy*.txt
    do
        tail -n +2 $fname >> results/{project_name}_classified.tsv
    done
   """.format(project_name=project_name)         
   
# FRIGGA will integrate the information of MOTU abundances and taxonomy assignment in a single table

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
        
input_files= []

input_files.append("results/{}_SWARM_seeds.fasta".format(project_name))
input_files.append("results/{}_All_MOTUs_classified.tsv".format(project_name))

output_files = []

output_files.append("results/{}_match_list.txt".format(project_name))
output_files.append("results/{}_Discarded_LULU.tsv".format(project_name))
output_files.append("results/{}_Curated_LULU.tsv".format(project_name))
output_files.append("results/{}_Deleted_LULU_fate.tsv".format(project_name))

gwf.target(
            name="loki_{}".format(project_name),
            inputs=input_files,
            outputs=output_files,
            cores=1,
            memory="32g",
            walltime="12:00:00",
        ) << """
            eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
            conda activate mjolnir
            cd results
            Rscript ../scripts/loki.r {project_name}
        """.format(project_name=project_name)         
                        
### Combine all the metadata files into one large file

input_files = glob('tmp/metadata*')
 
output_file = "results/metadata/{}_metadata.tsv".format(project_name)
    
gwf.target(
   name="metadata_{}".format(project_name),
   inputs=input_files,
   outputs=output_file,
   cores=1,
   memory="8g",
   walltime="00:10:00",
 ) << """
    head -n1 tmp/metadata_S111.tsv > results/metadata/COSQ_metadata.tsv
    for fname in tmp/metadata*
    do
        tail -n +2 $fname >> results/metadata/COSQ_metadata.tsv
    done
   """   

#Importantly, a few of the original sample names were incorrect with regard to the PSU replicate number.
#Also, the sample names from the original sequencing run and from resequencing were identical, as it was assumed they would be automatically merged. 
#However, MJOLNIR kept them separate, which is probably more transparent.
#Therefore, the metadata file was corrected with the script correct_metadata.R. See the file COSQ_metadata_reps.tsv for old and corrected sample names.

# RAGNAROC will change the names of the samples to recover the original names and will remove unnecessary columns

lulu = "results/{}_Curated_LULU.tsv".format(project_name)

classified = "results/{}_All_MOTUs_classified.tsv".format(project_name)

meta_file = "results/metadata/{}_metadata_new.tsv".format(project_name)

output_file = "results/{}_final_dataset.tsv".format(project_name)

gwf.target(
            name="ragnaroc_{}".format(project_name),
            inputs=[lulu, classified, meta_file],
            outputs=output_file,
            cores=1,
            memory="2g",
            walltime="00:10:00",
        ) << """
            eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
            conda activate mjolnir
            Rscript scripts/ragnaroc.r results/{project_name} {meta_file}
        """.format(project_name=project_name, meta_file=meta_file)        
