# apply filters
	# 1. same filter of minimum relative abundance applied to MOTUs
	# 2. remove sequences with different length than the modal for each MOTU, also only lengths of 313+-(n*3) are allowed
	# 3. remove numts, this is:
		# a. seqs with codon stops
		# b. for metazoan and 313bp length sequences, sequences with changes in 5 well preserved aa.
		# (the genetic code used among the mitochondrial is the one with less codon stops per read and for metazoans also the one with less aa changes in 5 well preserved positions)
	# 4. remove ESV with less than 5 reads

input_files = []
input_files.append("{}/ESV_Adcorr.csv".format(output_dir))
input_files.append("results/{}_metadata.tsv".format(project_name))
input_files.append("results/{}_Curated_LULU.tsv".format(project_name))

final_dir = "results/finalfiles"

output_files = "{}/XXX".format(final_dir) #what is the output?

gwf.target(
            name="filters_{}".format(project_name),
            inputs=input_files,
            outputs=output_files,
            cores=1,
            memory="196g",
            walltime="04:00:00",            
        ) << """
            eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
            conda activate mjolnir
            echo start filtering of ESV
        	Rscript scripts/Script2.1_DnoisE_ESV.R -i {output_dir}/ESV_Adcorr.csv -o ${final_dir} -a results/{project_name}_metadata.tsv -b results/{project_name}_Curated_LULU.tsv &
        	wait
	        echo finished
        """.format(project_name=project_name, output_dir=output_dir, final_dir=final_dir)