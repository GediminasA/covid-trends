tmp_dir = config["work_dir"]+"/tmp"
rez_dir = config["work_dir"]+"/rez"
lineages = config["lineages"].split()
reference_lineage = config["reference_lineage"]
juliaenv = "scripts/julia_modules/JuliaClusterAndTreeTools"

include: "snakefiles/initial_clustering.smk"
include: "snakefiles/s_protein_analysis.smk"
include: "snakefiles/mink_analysis.smk"
include: "snakefiles/downasmpling4timetree.smk"
include: "snakefiles/downasmpling4timetree_clustered.smk"
include: "snakefiles/summarising_analyses.smk"
include: "snakefiles/haplotypes_analysis.smk"
include: "snakefiles/binding_evaluation.smk"

rule setup_R:
    output:
        config["work_dir"]+"/RsetupDone.txt"
    conda:
        "envs/R_env.yaml"
    notebook:
        "notebooks/setupR.r.ipynb"

# install IJulia on main conda env.
rule setup_julia_environment_on_main_env:
    output:
        config["work_dir"]+"/JuliasetupDoneOnMain.txt"
    params:
        wdir = config["work_dir"]
    shell:
        """
            julia -e 'using Pkg; Pkg.instantiate()'
            julia -e 'using Pkg; Pkg.add("IJulia")'
            julia -e 'using Pkg; Pkg.add("Revise")'
            touch {params.wdir}/JuliasetupDoneOnMain.txt
        """

# installs IJulia on R conda env
rule setup_julia_environment:
    input:
        config["work_dir"]+"/RsetupDone.txt",
        config["work_dir"]+"/JuliasetupDoneOnMain.txt"
    output:
        config["work_dir"]+"/JuliasetupDone.txt"
    params:
        julia_project_location = "scripts/julia_modules/JuliaClusterAndTreeTools/",
	wdir = config["work_dir"]
    threads: 32
    conda:
        "envs/R_env.yaml"
    shell:
        """
            cwd=`pwd`
            cd {params.julia_project_location}
            julia -e 'using Pkg; Pkg.instantiate()'
            julia -e 'using Pkg; Pkg.add("IJulia")'
            julia -e 'using Pkg; Pkg.build("IJulia")'
            julia -e 'using Pkg; Pkg.add("RCall")'
            julia -e 'using Pkg; Pkg.build("RCall")'
            julia --project=. -e 'using Pkg; Pkg.instantiate()'
            julia --project=. -e 'using Pkg; Pkg.precompile()'
            julia --project=. -e 'using Pkg; Pkg.build()'
            cd $cwd
            touch {params.wdir}/JuliasetupDone.txt
        """

rule target:
    input:
        #rez = rez_dir+"/initial_set_minimap2_match_all.sam",
        report = rez_dir+"/pangolin_lineage_report.csv"
        #rez_dir+"/initial_set.fasta",

rule run_pangolin:
    input:
        fasta = config["sequences"]
    output:
        report = protected(rez_dir+"/pangolin_lineage_report.csv")
    params:
        out_dir = rez_dir,    
        tmp_dir = tmp_dir + "/pangolin_analysis", 
        alignment = "alignment_minimap.fasta",
        report = "pangolin_lineage_report.csv"
    conda:
        "envs/pangolin.yaml"
    threads: 80
    shell:
        """
            pangolin --add-assignment-cache
	    mkdir -p {params.tmp_dir}
            pangolin --use-assignment-cache -t {threads} --outfile {params.report} --outdir {params.out_dir}  --tempdir {params.tmp_dir}   --alignment-file {params.alignment}  --analysis-mode fast  --alignment {input.fasta}
        """ 

rule download_covid_dataset:
    output:
        directory(config["datasets_dir"]+"/sars-cov-2")
    conda:
        "envs/nextclade.yaml"
    shell:
        """
            nextclade dataset get --name 'sars-cov-2' --output-dir '{output[0]}'
        """

rule run_nextclade:
    input:
        fasta = config["sequences"],
        data = config["datasets_dir"]+"/sars-cov-2"
    output:
        alignment = protected(rez_dir+"/alignment_nextclasde.fasta.gz"),
        report = protected(rez_dir+"/nextclade_report.csv"),
        S_alignment = protected(rez_dir+"/gene_S.translation.fasta.gz") 
    log:
        "logs/nextclade.log"
    conda:
        "envs/nextclade.yaml"
    threads: 80
    params:
        out_dir = rez_dir 
    shell:
        """
       nextclade run --jobs {threads} --include-reference --output-translations \
       {params.out_dir}/gene_{{gene}}.translation.fasta.gz --output-csv\
       {output.report} --output-fasta \
       {output.alignment}\
       --input-dataset {input.data}\
       {input.fasta} &> {log}
        """
    
rule test:
    input:
        rez_dir+"/nextclade_report.csv"

rule rule:
    input:
        rez_dir+"/nextclade_report.csv"
