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
        aa = wdirb + "/setup_done.txt",
        test_data = "data/example_4binding.csv",
        data = rez_dir + "/select_mutants_against_antibodies.csv",
        structures = "data/covid-lt-contacts/pdb/antibodies/complexes"
    output:
        "binding_evaluator/analyzis_runs/rezults/promod_models_results_full.csv"
    params:
        wdir = "analyzis_runs"
    conda: "../binding_evaluator/envs/binding_evaluator.yaml"
    shell:
        """
            cd binding_evaluator
            snakemake --profile local --configfile configs/antibody.yaml --rerun-incomplete  \
            -k get_summary_of_binding \
            --config workdir={params.wdir}  
        """


    
    

rule run_binding:
    input:
        wdirb + "/setup_done.txt",
        wdirb + "/run_done.txt"
