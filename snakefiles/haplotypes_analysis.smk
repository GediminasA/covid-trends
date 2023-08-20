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
        haplotypes = work_dir + "/haplotype_info_initial_clean_Sprotein_analysis.csv" #nextclade outfile
    log:
        haplotypes = work_dir + "/haplotype_info_initial_clean_Sprotein_analysis.log.ipynb" #nextclade outfile
    params:
        environment = "scripts/julia_modules/JuliaClusterAndTreeTools",
        cwd =  os.getcwd()
    threads: 16
    conda:
        "../envs/R_env.yaml"
    shell:
        """export JULIA_NUM_THREADS={threads}
        papermill --prepare-only --progress-bar  notebooks/get_S_haplotypes_analysis.jl.ipynb {log} \
         -p wdir {params.cwd} -p outfile {output.haplotypes} \
         -p haplotype {input.haplotypes} \
         -p meta {input.meta} \
         -p contacts {input.contacts} \
         -p pangolin {input.pangolin} \
         -p environment {params.environment} """

