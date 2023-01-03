output_dir = "tmp/output_Ad_corr"

input_file = "{}/COSQ_000000015_Adcorr_denoised_ratio_d.csv".format(output_dir)

output_file = "{}/ESV_Adcorr.csv".format(output_dir)

gwf.target(
            name="csv_{}".format(project_name),
            inputs=input_file,
            outputs=output_file,
            cores=1,
            memory="16g",
            walltime="01:00:00",            
        ) << """
            head -n 1 {input_file} > {output_file}
    	    sed -i 's/id/motu,id/g' {output_file}
        """.format(project_name=project_name, input_file=infput_file, output_file=output_file)


output_dir = "tmp/output_Ad_corr"

cores = 6

for motu in MOTUS:
    
    input_files = ["{}/{}".format(motus_tab_dir,motu),"{}/{}.log".format(motus_tab_dir,motu)]
    
    output_file = "{}/{}_Adcorr_progress.txt".format(output_dir,motu)
        
    gwf.target(
                name=f"dnoise_{motu}",
                inputs=input_files,
                outputs=output_file,
                cores=6,
                memory="16g",
                walltime="4:00:00",
            ) << """
                awk 'NR>1' {output_dir}/{motu}_Adcorr_denoised_ratio_d.csv | awk -v var="$motu," '{print var $0;}' >>{output_dir}/ESV_Adcorr.csv
            """.format(motu=motu, output_dir=output_dir)

# apply filters
# 1. same filter of minimum relative abundance applied to MOTUs
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

lulu = "results/COSQ_Curated_LULU.tsv"

output_file = "{}/filter.log".format(output_dir)

gwf.target(
            name="filter_{}".format(project_name),
            inputs=[input_file,metadata,lulu]
            outputs=output_file,
            cores=1,
            memory="16g",
            walltime="01:00:00",            
        ) << """
            echo "start filtering of ESVs"
            Rscript t./scripts/Script2.1_DnoisE_ESV.R -i {input_file} -o {final_dir} -a {metadata} -b {lulu} &
            wait
	        echo "finished"
            echo "hello" > {output_dir}/filter.log
        """.format(project_name=project_name, input_file=input_file, final_dir=final_dir, metadata=metadata, lulu=lulu, output_file=output_file)
