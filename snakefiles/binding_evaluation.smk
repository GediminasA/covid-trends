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

rule run_antibodies_evaluation:
    input:
        setup_marker = wdirb + "/setup_done.txt",
        test_data = "data/example_4binding.csv",
        data = rez_dir + "/select_mutants_against_antibodies.csv",
        structures = "data/covid-lt-contacts/pdb/antibodies/complexes"
    output:
        "binding_evaluator/analysis_runs/rezults/promod_models_results_full.csv"
    log:
         wdirb + "/binding_evaluation.log",
    params:
        wdir = "analysis_runs",
        starting_dir = os.getcwd()  
    conda: "../binding_evaluator/envs/binding_evaluator.yaml"
    shell:
        """
            cd binding_evaluator
            script -efq -c '\
            snakemake --profile local --configfile configs/antibody.yaml --rerun-incomplete  \
            -k get_summary_of_binding \
            --config workdir={params.wdir} \
            '  |& tee -i {params.starting_dir}/{log} 
            echo "LOG OF BINDING EVALUATION'S WORKFLOW: {params.starting_dir}/{log}"
        """


    
    

rule run_binding:
    input:
        wdirb + "/setup_done.txt",
        wdirb + "/run_done.txt"
