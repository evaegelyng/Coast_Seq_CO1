from gwf import Workflow

project_name = "COSQ"

gwf = Workflow(defaults={"account": "edna"}) 

# apply filters
# 1. same filter of minimum relative abundance applied to MOTUs. EES: This step was omitted, as a similar cleaning was already done before denoising
# 2. remove sequences with different length than the modal for each MOTU, also only lengths of 313+-(n*3) are allowed
# 3. remove numts, this is:
	# a. seqs with codon stops
	# b. for metazoan and 313bp length sequences, sequences with changes in 5 well preserved aa.
	# (the genetic code used among the mitochondrial is the one with less codon stops per read and for metazoans also the one with less aa changes in 5 well preserved positions)
# 4. remove ESV with less than 5 reads

output_dir = "tmp/output_Ad_corr"

final_dir = "tmp/finalfiles"

input_file = "{}/ESV_Adcorr.csv".format(output_dir)

metadata = "results/{}_metadata.tsv".format(project_name)

lulu = "results/{}_Curated_LULU.tsv".format(project_name)

output_file = "{}/filter.log".format(output_dir)

gwf.target(
            name="filter_{}".format(project_name),
            inputs=[input_file,metadata,lulu],
            outputs=output_file,
            cores=1,
            memory="16g",
            walltime="01:00:00",            
        ) << """
            echo "start filtering of ESVs"
            Rscript ../scripts/Script2.1_DnoisE_ESV.R -i {input_file} -o {final_dir} -a {metadata} -b {lulu} &
            wait
	        echo "finished"
            echo "hello" > {output_dir}/filter.log
        """.format(project_name=project_name, input_file=input_file, final_dir=final_dir, metadata=metadata, lulu=lulu, output_file=output_file, output_dir=output_dir)
