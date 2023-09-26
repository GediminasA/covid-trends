wdirb = tmp_dir+"/binding" # work dir for binding evals

rule setup_binding_submodule:
    output:
        wdirb + "/setup_done.txt"
    log:
        wdirb + "/setup.log"
    conda: "../binding_evaluator/envs/binding_evaluator.yaml"
    shell:
        """
            wd=`pwd`
            cd binding_evaluator
            ls
            snakemake --profile local --configfile configs/antibody.yaml  -k get_summary_of_binding
            cd $wd
            touch {output}
        """

rule get_config_4_antibodies_evaluation:
    input:
        setup_marker = wdirb + "/setup_done.txt",
        test_data = "data/example_4binding.csv",
        data = rez_dir + "/select_mutants_against_antibodies.csv",
        structures = "data/covid-lt-contacts/pdb/antibodies/complexes"
    output:
        "binding_evaluator/configs/antibody_runs.yaml"
    params:
        wdir = "analysis_of_antibodies_runs",
        starting_dir = os.getcwd()  
    run:
        confcont = f"""structures_folder: data/antibody 
preprocessed_structures: {params.starting_dir}/{input.structures} #data/antibody_test_structures 
pdb_stems: 7LQV
mutants: {params.starting_dir}/{input.data} #data/example_mutations.csv
mutants_templates: data/example_mations_templates.fasta 
workdir: {params.wdir}
conformers_to_evaluate: 0 """
        with ( open(output[0], "w") as conf ):
            conf.write(confcont)

        

rule run_antibodies_evaluation:
    input:
        setup_marker = wdirb + "/setup_done.txt",
        conf = "binding_evaluator/configs/antibody_runs.yaml"
    output:
        run_scripts = "binding_evaluator/analysis_of_antibodies_runs.sh",
        #rezults = "binding_evaluator/analysis_of_antibodies_runs/rezults/promod_models_results_full.csv"
    log:
         wdirb + "/binding_evaluation.log",
    params:
        wdir = "analysis_of_antibodies_runs",
        starting_dir = os.getcwd(), 
        localconf = "configs/antibody_runs.yaml",
    conda: "../binding_evaluator/envs/binding_evaluator.yaml"
    shell:
        """
            cd binding_evaluator
            echo snakemake --profile local \
            --rerun-incomplete --configfile {params.localconf} -k -f get_summary_of_binding \
            > {params.starting_dir}/{output.run_scripts} 
            script=`basename {output.run_scripts}`
            chmod a+x $script
            echo running $script in directory `pwd`
            #script -efq -c " ./$script"  |& tee -i {params.starting_dir}/{log} 
            echo "LOG OF BINDING EVALUATION'S WORKFLOW: {params.starting_dir}/{log}"
        """


    
    

rule get_config_4_ace2_evaluation:
    input:
        setup_marker = wdirb + "/setup_done.txt",
        test_data = "data/example_4binding.csv",
        data = rez_dir + "/select_mutants_against_ace2.csv",
        structures = "data/covid-lt-contacts/pdb/ACE2/complexes"
    output:
        "binding_evaluator/configs/ace2_runs.yaml"
    params:
        wdir = "analysis_of_ace2_runs",
        starting_dir = os.getcwd()  
    run:
        confcont = f"""structures_folder: data/ace2 
preprocessed_structures: {params.starting_dir}/{input.structures} #data/antibody_test_structures 
pdb_stems: 7LQV
mutants: {params.starting_dir}/{input.data} #data/example_mutations.csv
mutants_templates: data/example_mations_templates.fasta 
workdir: {params.wdir}
conformers_to_evaluate: 0 """
        with ( open(output[0], "w") as conf ):
            conf.write(confcont)

rule run_ace2_evaluation:
    input:
        setup_marker = wdirb + "/setup_done.txt",
        conf = "binding_evaluator/configs/ace2_runs.yaml"
    output:
        run_scripts = "binding_evaluator/analysis_of_ace2_runs.sh",
        #rezults = "binding_evaluator/analysis_of_antibodies_runs/rezults/promod_models_results_full.csv"
    log:
         wdirb + "/binding_evaluation.log",
    params:
        wdir = "analysis_of_ace2_runs",
        starting_dir = os.getcwd(), 
        localconf = "configs/ace2_runs.yaml",
    conda: "../binding_evaluator/envs/binding_evaluator.yaml"
    shell:
        """
            cd binding_evaluator
            echo snakemake --profile local \
            --rerun-incomplete --configfile {params.localconf} -k -f get_summary_of_binding \
            > {params.starting_dir}/{output.run_scripts} 
            script=`basename {output.run_scripts}`
            chmod a+x $script
            echo running $script in directory `pwd`
            #script -efq -c " ./$script"  |& tee -i {params.starting_dir}/{log} 
            echo "LOG OF BINDING EVALUATION'S WORKFLOW: {params.starting_dir}/{log}"
        """

rule run_binding:
    input:
        #"binding_evaluator/analysis_of_antibodies_runs/rezults/promod_models_results_full.csv"
       # "binding_evaluator/configs/antibody_run.yaml"
        "binding_evaluator/analysis_of_antibodies_runs.sh",
        "binding_evaluator/analysis_of_ace2_runs.sh",
        #wdirb + "/setup_done.txt",
        #wdirb + "/run_done.txt"
