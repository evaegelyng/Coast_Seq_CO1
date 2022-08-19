#ESV Addcor

with open(MOTUS2RUN, 'r') as fp:
    read = fp.readlines() 
    lines = len(read) # Result: 216269, covering read[0] to reads[216268]. Seems like the
    #result from bash was incorrect by one line?

for i in range(1,len(read)):
    motu = read[i].strip()
    
    input_files = []
    input_files.append("{}/{}_Adcorr_denoised_ratio_d.csv".format(output_dir,motu))
    input_files.append(MOTUS2RUN)

    output_file = "{}/esv_done_{}.txt".format(output_dir,i) # Only added the first MOTU tab file as an output
        
    gwf.target(
                name="esv_{}_{}".format(project_name, i),
                inputs=input_files,
                outputs=output_file,
                cores=1,
                memory="4g",
                walltime="4:00:00",
            ) << """
                scripts/esv.sh {i} {MOTUS2RUN} {output_dir}
                echo ${i} > {output_dir}/esv_done_{i}.txt
            """.format(i=i, MOTUS2RUN=MOTUS2RUN, output_dir=output_dir)