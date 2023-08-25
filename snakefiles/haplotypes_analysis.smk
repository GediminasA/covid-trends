work_dir = config["work_dir"]+"/tmp/haplotypes"



rule target_hap:
    input:
        work_dir + "/metadata_info_initial_clean.csv",
        work_dir + "/haplotype_info_initial_clean.csv",
        work_dir + "/haplotype_info_initial_clean_Sprotein.csv.gz"

rule get_metadata_pango_clean:
    input:
        pangolin = rez_dir+"/pangolin_lineage_report.csv",
        metadata = config["meta"]
    output:
        cleanmeta = work_dir + "/metadata_info_initial_clean.csv"
    log:
        notebook = work_dir + "/metadata_info.r.ipynb.log"
    conda:
        "../envs/R_env.yaml"
    notebook:
        "../notebooks/get_metadata_pango_clean.r.ipynb"

rule clean_haplotypes_initial:
    input:
        nextclade =  rez_dir+"/nextclade_report.csv",
        cleanmeta = work_dir + "/metadata_info_initial_clean.csv"
    output:
        haplotypes = work_dir + "/haplotype_info_initial_clean.csv"
    log:
        haplotypes = work_dir + "/haplotype_info_initial_clean.log"
    conda:
        "../envs/R_env.yaml"
    notebook:
        "../notebooks/get_haplotypes_initial_clean.r.ipynb"


rule choose_S_protein_haplotypes:
    input:
        haplotypes = work_dir + "/haplotype_info_initial_clean.csv"
    output:
        haplotypes = work_dir + "/haplotype_info_initial_clean_Sprotein.csv.gz"
    log:
        haplotypes = work_dir + "/haplotype_info_initial_clean_Sprotein.log.ipynb" #nextclade outfile
    params:
        environment = "scripts/julia_modules/JuliaClusterAndTreeTools",
        cwd =  os.getcwd()
    shell:
        "papermill --progress-bar  notebooks/get_S_haplotypes.jl.ipynb {log} -p wdir {params.cwd} -p nextclade {input} -p outfile {output} -p environment {params.environment}"
        #"papermill   --prepare-only -- {params.cwd} --progress-bar  notebooks/get_S_haplotypes.jl.ipynb {log} -p wdir {params.cwd} -p nextclade {input} -p outfile {output} -p environment {params.environment}"

rule review_contact_data:
    input:
        contacts = config["antibody_contacts_data"]
    output:
        rez_dir + "/contacts_review.svg",
        rez_dir + "/contacts_per_cluster.csv"
    conda:
        "../envs/R_env.yaml"
    notebook:
        "../notebooks/review_contacts.r.ipynb"

rule analyse_S_protein_haplotype:
    input:
        haplotypes = work_dir + "/haplotype_info_initial_clean_Sprotein.csv.gz",
        meta = work_dir + "/metadata_info_initial_clean.csv",
        contacts = rez_dir + "/contacts_per_cluster.csv",
        pangolin = rez_dir+"/pangolin_lineage_report.csv",
    output:
        review = rez_dir + "/info_on_cluster_and_contacts_review.png",
        all_trend = rez_dir + "/info_on_cluster_and_contacts_over_time_all.svg",
        lt_trend = rez_dir + "/info_on_cluster_and_contacts_over_time_lt.svg",
        taus = rez_dir + "/info_on_cluster_and_contacts_over_time.csv",
        haplotypes = work_dir + "/haplotype_info_lt_per_month.csv"
    log:
        haplotypes = work_dir + "/haplotype_info_lt_per_month.log.ipynb"
    params:
        environment = "scripts/julia_modules/JuliaClusterAndTreeTools",
        cwd =  os.getcwd()
    threads: 16
    conda:
        "../envs/R_env.yaml"
    shell:
        """export JULIA_NUM_THREADS={threads}
        papermill  --progress-bar  notebooks/get_S_haplotypes_analysis.jl.ipynb {log} \
         -p wdir {params.cwd} -p outfile {output.haplotypes} \
         -p haplotype {input.haplotypes} \
         -p meta {input.meta} \
         -p contacts {input.contacts} \
         -p pangolin {input.pangolin} \
         -p environment {params.environment} \
         -p review {output.review} \
         -p all_trend {output.all_trend} \
         -p lt_trend {output.lt_trend} \
         -p taus {output.taus} \
         -p haplotypes {output.haplotypes}
         """

rule select_4_analysis:
    input:
        df_file = work_dir + "/haplotype_info_lt_per_month.csv",
        haplotype = work_dir + "/haplotype_info_initial_clean_Sprotein.csv.gz",
        contacts = rez_dir + "/contacts_per_cluster.csv",
        antibody_contacts_data_pdbwise = config["antibody_contacts_data_pdbwise"],
        pdb_cluster_assignments = config["antibody_contacts_data_cluster_assignments"],
    output:
        additional_data = rez_dir + "/select_mutants_against_antibodies_additional_info.csv",
        chosen_data = rez_dir + "/select_mutants_against_antibodies.csv",
    params:
        max_sructures_to_test_per_month4antib = config["max_sructures_to_test_per_month4antib"]
    log:
        notebook = rez_dir + "/select_mutants_against_antibodies.ipynb.log",
    conda:
        "../envs/R_env.yaml"
    notebook:
        "../notebooks/select_mutants_against_antibodies.r.ipynb"

rule run_haplotypes_analysis:
    input:
        rez_dir + "/select_mutants_against_antibodies.csv",
