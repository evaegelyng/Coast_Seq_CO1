# Run DnoisE. Remember to input entropy values to dnoise.sh

output_dir = "tmp/output_Ad_corr"

with open(MOTUS2RUN, 'r') as fp:
    read = fp.readlines() 
    lines = len(read) # Result: 216269, covering read[0] to reads[216268]. Seems like the
    #result from bash was incorrect by one line?

for i in range(1,len(read)):
    #input_file = "{}/{}_{}.csv".format(motus_tab_dir,project_name,var)
    input_file = MOTUS2RUN
    output_file = "{}/done_{}.txt".format(output_dir,i) # Only added the first MOTU tab file as an output
        
    gwf.target(
                name="dnoise_{}_{}".format(project_name, i),
                inputs=input_file,
                outputs=output_file,
                cores=6,
                memory="16g",
                walltime="4:00:00",
            ) << """
                mkdir -p {output_dir}
                scripts/dnoise.sh {i} {MOTUS2RUN} {output_dir} {DnoisE_dir} {motus_tab_dir} {cores}
                echo ${i} > {output_dir}/done_{i}.txt
            """.format(i=i, MOTUS2RUN=MOTUS2RUN, output_dir=output_dir, DnoisE_dir=DnoisE_dir, motus_tab_dir=motus_tab_dir, cores=cores)
