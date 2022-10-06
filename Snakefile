tmp_dir = config["work_dir"]+"/tmp"
rez_dir = config["work_dir"]+"/rez"


rule target:
    input:
        rez = rez_dir+"/initial_set_minimap2_match_all.sam",
        report = rez_dir+"/pangolin_lineage_report.csv"
        #rez_dir+"/initial_set.fasta",

rule get_subset_of_animals:
    output:
        rez_dir+"/initial_set_ids.txt"
    params:
        meta = config['nextclade_data_lt_mink']
    notebook:
        "notebooks/get_initial_set_ids.r.ipynb"

rule get_fasta:
    input:
        ids = "audines0520/rez/initial_set_ids.txt",
        seqs = config["sequences_lt_mink"]
    output:
        rez_dir+"/initial_set.fasta"
    shell:
        """
            seqkit grep -n -f {input.ids} -w 0 {input.seqs} > {output}
        """

rule get_hisat_index:
    input: "{stem}.fasta"
    output: "{stem}_hisat2.1.ht2"
    params: stem = "{stem}_hisat2"
    threads: 80
    conda:
        "envs/analysis.yaml"
    shell: "hisat2-build -p {threads} {input} {params.stem}"    

rule remove_n_subs:
    input:
        reference = rez_dir+"/initial_set.fasta",
    output:
        reference = rez_dir+"/initial_set_woN.fasta",
    conda:
        "envs/analysis.yaml"
    shell:
       " seqkit seq -g -G N -i {input} -o {output} "

rule remove_n_alls:
    input:
        fasta = config["sequences"]
    output:
        raw_woN = temp(rez_dir+"/tmp/lineages_all_woN.fasta")
    conda:
        "envs/analysis.yaml"
    shell:
       " seqkit seq -g -G N -i {input} -o {output} "


rule find_closestv2:
    input:
        reference = rez_dir+"/initial_set_woN.fasta",
        reference_index = rez_dir+"/initial_set_hisat2.1.ht2",
        fasta = temp(rez_dir+"/tmp/lineages_all_woN.fasta")
    output:
        rez = rez_dir+"/initial_set_minimap2_match_all.sam"
    params:
        index = rez_dir+"/initial_set_hisat2"
    threads: 80
    shell:
        """
        hisat2 -p {threads}  -x {params.index}  -f {input.fasta}   --no-spliced-alignment  -k 2 > {output}
        """
rule run_pangolin:
    input:
        fasta = config["sequences"]
    output:
        alignment = rez_dir+"/alignment_minimap.fasta",
        report = rez_dir+"/pangolin_lineage_report.csv"
    params:
        calculation_dir = rez_dir + "pangolin_analysis",    
        tmp_dir = rez_dir + /"pangolin_analysis" 
    conda:
        "envs/pangolin.yaml"
    threads: 80
    shell:
        """
          pangolin -t {threads} --outfile {output.report} --outdir {params.calculation_dir}  --tempdir {params.tmp_dir}   --alignment-file {output.alignment}  --analysis-mode fast  --alignment {input.fasta}
        """ 
 

