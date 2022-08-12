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
