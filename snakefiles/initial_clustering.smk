rule clusteringI:    
    input:
        #expand(rez_dir+ "/lineages/{id}/ids.txt",id=lineages)
        #expand(rez_dir + "/lineages/{id}/alignment_nextclade_woG.fasta.gz",id=lineages),
        #expand(rez_dir + "/lineages/{id}/alignment_nextclade_woG_derep1.fasta.uc",id=lineages),
        #expand(rez_dir + "/lineages/{id}/alignment_nextclade_woG_derep1_mappings.ids",id=lineages),
        #expand(rez_dir + "/lineages/{id}/alignment_nextclade_derep.pssm",id=lineages)
        #rez_dir + "/lineages/alignment_nextclade_refseq_shortid.fasta.gz",
        expand(rez_dir + "/lineages/{id}/alignment_nextclade_filtered_0.8_20_woG_derep1_swarm.fasta",id=lineages)
        #expand(rez_dir + "/lineages/{id}/alignment_nextclade_filtered_0.8_20_woG_derep1.fasta.gz",id=lineages)
rule get_ids:
    input:
        pangolin = rez_dir+"/pangolin_lineage_report.csv",
        nextclade = rez_dir+"/nextclade_report.csv",
        meta = config["meta"]
    output:
        ids = rez_dir + "/lineages/{id}/ids.txt",
        data = rez_dir + "/lineages/{id}/data.csv",
    params:
        id = "{id}"
    threads: 12
    notebook:
        "notebooks/filter_out_lioneage.r.ipynb"

rule extract_alignment:
    input:
        ids = rez_dir + "/lineages/{id}/ids.txt",
        alignment = rez_dir+"/alignment_nextclasde.fasta.gz",
        refclean = rez_dir + "/lineages/alignment_nextclade_refseq_shortid.fasta.gz",
    output:
        aln = rez_dir + "/lineages/{id}/alignment_nextclade.fasta.gz",
    threads: 12
    shell:
        """
            
            seqkit grep -w 0 -n -f {input.ids} {input.alignment} -o {output.aln}
            cat {input.refclean} >>   {output.aln}
        """

rule extract_reference_sequence:
    input:
        alignment = rez_dir+"/alignment_nextclasde.fasta.gz",
    output:
        ref = rez_dir + "/lineages/alignment_nextclade_refseq.fasta.gz",
    params:
        id = config["reference_id"]
    threads: 12
    shell:
        "seqkit grep -w 0 -n -p '{params.id}' {input.alignment} -o {output.ref}  "


rule clean_name_4_reference_sequence:
    input:
        ref = rez_dir + "/lineages/alignment_nextclade_refseq.fasta.gz",
    output:
        refclean = rez_dir + "/lineages/alignment_nextclade_refseq_shortid.fasta.gz",
    shell:
        "seqkit seq -w 0 -i {input} -o {output}  "

# run filter indudually leaving only sequences without Ns
rule individual_filter:
    input:
        in_file_name = rez_dir + "/lineages/{id}/alignment_nextclade.fasta.gz",
    output:
        output_file_name = rez_dir + "/lineages/{id}/alignment_nextclade_filtered_{k}_{t,\d+}.fasta.gz", 
        output_file_nonconspoz_name = rez_dir + "/lineages/{id}/alignment_nextclade_filtere_nonconspoz_{k}_{t,\d+}.fasta.gz",
        rows_after_filter = rez_dir + "/lineages/{id}/alignment_nextclade_filtered_keptrows_{k}_{t}.txt", 
        columns_after_filter = rez_dir + "/lineages/{id}/alignment_nextclade_filtered_keptcolumns_{k}_{t}.txt", 
        conservation_data =  rez_dir + "/lineages/{id}/alignment_nextclade_filtered_consData_{k}_{t}.csv", 
    params:
        keepPositionsFraction = "{k}",#0.8, # 0.7 0.8 0.9 0.95
        topNhaving_sequences_to_remove ="{t}", # 10, # 10 20 10 5 2
        reference_seq_id = "MN908947",
        country = "Lithuania",
    threads: 1000
    log: rez_dir + "/lineages/{id}/alignment_nextclade_filtered_parsed_{k}_{t}.ipynb"
    shell:
        """
        papermill scripts/julia_modules/JuliaClusterAndTreeTools/notebooks/test_aln.ipynb \
        {log} \
        -p in_file_name {input.in_file_name} \
        -p  output_file_name {output.output_file_name} \
        -p output_file_nonconspoz_name {output.output_file_nonconspoz_name} \
        -p rows_after_filter {output.rows_after_filter} \
        -p columns_after_filter {output.columns_after_filter} \
        -p conservation_data {output.conservation_data} \
        -p keepPositionsFraction {params.keepPositionsFraction} \
        -p topNhaving_sequences_to_remove {params.topNhaving_sequences_to_remove} \
        -p reference_seq_id {params.reference_seq_id} \
        -p country {params.country} 
        """

rule remove_s_with_N:
    input:
        "{stem}.fasta.gz"
    output:
        "{stem}_woN.fasta.gz"
    threads: 12
    shell:
        "bbduk.sh t={threads} ordered=t in={input} out={output} maxns=0  "

rule remove_gaps:
    input:
        "{stem}.fasta.gz"
    output:
        "{stem}_woG.fasta.gz"
    threads: 12
    shell:
        """
            seqkit seq -w 0 -g  {input}  | gzip --stdout >  {output}
        """




rule rmident1:
    input:
        "{stem}.fasta.gz",
    output:
        fa = "{stem}_derep1.fasta",
        uc = "{stem}_derep1.fasta.uc"
    params:
    conda: "../envs/clustering_tools.yaml"
    shell:
        '''
        vsearch    --derep_fulllength   {input}    --sizeout   --fasta_width 0  --output {output[0]} --uc {output.uc}
        '''

rule swarm:
    input:
        "{stem}.fasta",
    output:
        fa = "{stem}_swarm.fasta",
        uc = "{stem}_swarm.fasta.swarminfo"
    params:
    conda: "../envs/clustering_tools.yaml"
    threads: 1000
    shell:
        '''
        swarm -f --threads {threads} -z   -w {output[0]} -r -o {output.uc} {input}
        '''

rule parse_uc:
    input:
        uc = "{stem}.fasta.uc",
        setupmark = config["work_dir"]+"/RsetupDone.txt"
    output:
        data = "{stem}_mappings.rds",
        rep_ids = "{stem}_mappings.ids",
    conda:
        "../envs/R_env.yaml"
    notebook:
        "../notebooks/parse_uc.r.ipynb"


rule extract_alignment_derep:
    input:
        ids = rez_dir + "/lineages/{id}/alignment_nextclade_woG_derep1_mappings.ids",
        alignment = rez_dir + "/lineages/{id}/alignment_nextclade.fasta.gz",
    output:
        aln = rez_dir + "/lineages/{id}/alignment_nextclade_derep.fasta.gz",
    threads: 12
    shell:
        "seqkit grep -w 0 -n -f {input.ids} {input.alignment} -o {output.aln}  "

rule get_pssm:
    input:    
        rez_dir + "/lineages/{id}/alignment_nextclade_derep.fasta.gz"
    output:
        rez_dir + "/lineages/{id}/alignment_nextclade_derep.pssm"
    conda: "../envs/clustering_tools.yaml"
    threads: 12
    shell:
        """
            goalign compute pssm -i {input} --auto-detect  -n 0 -t {threads} > {output}
        """
