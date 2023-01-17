cleanup_keep_fraction_of_columns = str(config["cleanup_keep_fraction_of_columns"])
cleanup_remove_proc_mostN_having = str(config["cleanup_remove_proc_mostN_having"])
reference_lineage = config["reference_lineage"]
lineages4summary = set(lineages) - set([config["reference_lineage"]])
rule clusteringI:    
    input:
        #expand(rez_dir+ "/lineages/{id}/ids.txt",id=lineages)
        #expand(rez_dir + "/lineages/{id}/alignment_nextclade_woG.fasta.gz",id=lineages),
        #expand(rez_dir + "/lineages/{id}/alignment_nextclade_woG_derep1.fasta.uc",id=lineages),
        #expand(rez_dir + "/lineages/{id}/alignment_nextclade_woG_derep1_mappings.ids",id=lineages),
        #expand(rez_dir + "/lineages/{id}/alignment_nextclade_derep.pssm",id=lineages)
        #rez_dir + "/lineages/alignment_nextclade_refseq_shortid.fasta.gz",
        expand(rez_dir + "/lineages/{id}/alignment_nextclade_filtered_{fr}_{topn}_woG_derep1.fasta",
        id=lineages,
        fr= cleanup_keep_fraction_of_columns,
        topn=cleanup_remove_proc_mostN_having)
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

rule get_ids_lt:
    input:
        data = rez_dir + "/lineages/{id}/data.csv",
    output:
        ids = rez_dir + "/lineages/{id}/ids_lt.txt",
    params:
        id = "{id}"
    threads: 12
    conda:
        "../envs/R_env.yaml"
    notebook:
        "notebooks/filter_out_ltseqs.r.ipynb"


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
        papermill  scripts/julia_modules/JuliaClusterAndTreeTools/notebooks/test_aln.ipynb \
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

rule swarmF:
    input:
        "{stem}.fasta",
    output:
        fa = "{stem}_swarmF.fasta",
        uc = "{stem}_swarmF.fasta.swarminfo"
    params:
        gap_opening_penalty = 100,
        gap_extension_penalty = 100
    conda: "../envs/clustering_tools.yaml"
    threads: 1000
    shell:
        '''
        swarm -f --threads {threads} -z   -w {output[0]} -r -o {output.uc} {input} -e {params.gap_extension_penalty} -g {params.gap_opening_penalty}
        '''

rule swarm:
    input:
        "{stem}.fasta",
    output:
        fa = "{stem}_swarm.fasta",
        uc = "{stem}_swarm.fasta.swarminfo",
        stru = "{stem}_swarm.fasta.internstr"
    params:
    conda: "../envs/clustering_tools.yaml"
    threads: 1000
    shell:
        '''
        swarm -f --threads {threads} -z -d 1   -w {output[0]} -r -i {output.stru} -o {output.uc} {input}
        '''

# rule swarm:
#     input:
#         "{stem}.fasta",
#     output:
#         fa = "{stem}_swarm.fasta",
#         uc = "{stem}_swarm.fasta.swarminfo"
#     params:
#         gap_opening_penalty = 100,
#         gap_extension_penalty = 100
#     conda: "../envs/clustering_tools.yaml"
#     threads: 1000
#     shell:
#         '''
#         swarm -f --threads {threads} -z   -w {output[0]} -r -o {output.uc} {input} -e {params.gap_extension_penalty} -g {params.gap_opening_penalty}
#         '''
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

# pairwise analysis
rule get_kept_columns:
    input:
        ref_kept = rez_dir + "/lineages/"+reference_lineage+"/alignment_nextclade_filtered_keptcolumns_"+cleanup_keep_fraction_of_columns+"_"+cleanup_remove_proc_mostN_having+".txt",
        id_kept = rez_dir + "/lineages/{id}/alignment_nextclade_filtered_keptcolumns_"+cleanup_keep_fraction_of_columns+"_"+cleanup_remove_proc_mostN_having+".txt",
    params:
        id = "{id}",
        ref = reference_lineage
    output:
        id_kept = rez_dir + "/lineages/{id}/kept_columns_reference_based.txt",
    log:
        notebook = rez_dir + "/lineages/{id}/kept_columns_reference_based.r.ipynb.log",
    conda:
        "../envs/R_env.yaml"
    notebook:
        "../notebooks/get_common_rows.r.ipynb"

        #ref_kept = rez_dir + "/lineages/{id}/alignment_nextclade_filtered_keptrows_{cleanup_keep_fraction_of_columns}_{cleanup_remove_proc_mostN_having}.txt",

rule get_subset_alignments:
    input:
        columns_chosen = rez_dir + "/lineages/{id}/kept_columns_reference_based.txt",
        rows_chosen_id = rez_dir + "/lineages/{id}/alignment_nextclade_filtered_keptrows_"+cleanup_keep_fraction_of_columns+"_"+cleanup_remove_proc_mostN_having+".txt",
        rows_chosen_reference = rez_dir + "/lineages/"+reference_lineage+"/alignment_nextclade_filtered_keptrows_"+cleanup_keep_fraction_of_columns+"_"+cleanup_remove_proc_mostN_having+".txt",
        aln_id = rez_dir + "/lineages/{id}/alignment_nextclade.fasta.gz",
        aln_ref = rez_dir + "/lineages/"+reference_lineage+"/alignment_nextclade.fasta.gz",
    output:
        output_id = rez_dir + "/lineages/{id}/pair_id.fasta.gz",
        output_ref = rez_dir + "/lineages/{id}/pair_ref.fasta.gz"
    threads: 16
    log: rez_dir + "/lineages/{id}/pair_extract_names_columns_aln.ipynb"
    shell:
        """
        papermill scripts/julia_modules/JuliaClusterAndTreeTools/notebooks/extract_names_columns_aln.ipynb \
        {log} \
        -p columns_chosen {input.columns_chosen} \
        -p rows_chosen_id {input.rows_chosen_id} \
        -p rows_chosen_reference {input.rows_chosen_reference} \
        -p aln_id {input.aln_id} \
        -p aln_ref {input.aln_ref} \
        -p output_id {output.output_id} \
        -p output_ref {output.output_ref} 
        """


rule get_subset_alignments_common_cleaned:
    input:
        columns_chosen = expand(rez_dir + "/lineages/{id}/kept_columns_reference_based.txt", id = lineages),
        rows_chosen_id = rez_dir + "/lineages/{id}/alignment_nextclade_filtered_keptrows_"+cleanup_keep_fraction_of_columns+"_"+cleanup_remove_proc_mostN_having+".txt",
        aln_id = rez_dir + "/lineages/{id}/alignment_nextclade.fasta.gz",
    output:
        output_id = rez_dir + "/lineages/{id}/common_id.fasta.gz",
    threads: 1000
    log: rez_dir + "/lineages/{id}/pair_extract_names_columns_aln_common.ipynb"
    shell:
        """
        papermill scripts/julia_modules/JuliaClusterAndTreeTools/notebooks/pair_extract_names_columns_aln_common.ipynb \
        {log} \
        -p columns_chosen "{input.columns_chosen}" \
        -p rows_chosen_id "{input.rows_chosen_id}" \
        -p aln_id {input.aln_id} \
        -p output_id {output.output_id} 
        """
#replaces gaps with 'a'
rule replace_gaps:
    input:
        "{stem}.fasta.gz"
    output:
        "{stem}_Gap2a.fasta.gz"
    log:
        "{stem}_Gap2a.log"
    shell:
        """
        papermill scripts/julia_modules/JuliaClusterAndTreeTools/notebooks/replace_gaps.ipynb \
        {log} \
        -p in_file_name {input} \
        -p output_file_name {output} 
        """

#analyse clusterings
rule amalyze_pair_derep:
    input:
        data = rez_dir + "/lineages/{id}/data.csv",
        data_ref = rez_dir + "/lineages/"+reference_lineage+"/data.csv",
        id = rez_dir + "/lineages/{id}/pair_id_Gap2a_derep1.fasta.uc",
        ref = rez_dir + "/lineages/{id}/pair_ref_Gap2a_derep1.fasta.uc",
        setupmark = config["work_dir"]+"/RsetupDone.txt"
    log:
        notebook = rez_dir + "/lineages/{id}/pair_ref_derep_data.r.ipynb"
    output:
        ref = rez_dir + "/lineages/{id}/pair_ref_derep_data.csv"
    conda:
        "../envs/R_env.yaml"
    notebook:
        "../notebooks/analyse_derep_pair.r.ipynb"

rule amalyze_common_derep:
    input:
        datas = expand(rez_dir + "/lineages/{id}/data.csv", id = lineages),
        ids = expand(rez_dir + "/lineages/{id}/common_id_derep1.fasta.uc", id = lineages),
        setupmark = config["work_dir"]+"/RsetupDone.txt"
    log:
        notebook = rez_dir + "/lineages/common/common_ref_derep_data.r.ipynb"
    output:
        ref = rez_dir + "/lineages/common/common_ref_derep_data.csv"
    conda:
        "../envs/R_env.yaml"
    notebook:
        "../notebooks/analyse_derep_common.r.ipynb"

rule sumap_pair_derep:
    input:
        refs= expand(rez_dir + "/lineages/{id}/pair_ref_derep_data.csv", id = lineages4summary)
    output:
        rez_dir + "/pairwise_derep_results.csv" 
    log:
        notebook = rez_dir + "/pairwise_derep_results.csv"
    conda:
        "../envs/R_env.yaml"
    notebook:
        "../notebooks/summarise_derep_pair.r.ipynb"

     
rule cleanup_swarm_clusters:
    input:
        id = rez_dir + "/lineages/{id}/pair_id_Gap2a_derep1_swarm.fasta.swarminfo",
    output:
        id = rez_dir + "/lineages/{id}/pair_id_Gap2a_derep1_swarm.fasta.swarminfo.cleaned.txt"
    conda:
        "../envs/R_env.yaml"
    notebook:
        "../notebooks/clean_swarm.r.ipynb"

rule split_swarm_clusters:
    input:
        swarm =  rez_dir + "/lineages/{id}/pair_id_Gap2a_derep1_swarm.fasta.swarminfo.cleaned.txt",
        fasta = rez_dir + "/lineages/{id}/pair_id_Gap2a_derep1.fasta",
    output:
        swarm_single =  rez_dir + "/lineages/{id}/pair_id_Gap2a_derep1_swarm.fasta.swarminfo.cleaned_single.txt",
    run:
        print(input[0])
        print(input[1])
        print(output[0])

        




rule test2:
    input:
        #expand(rez_dir + "/lineages/{id}/kept_rows.txt",id=["BA.2"])
        expand(rez_dir + "/lineages/{id}/pair_id_Gap2a_derep1_swarm.fasta",id=["Q.1","B.1.1.7"]),
        expand(rez_dir + "/lineages/{id}/pair_ref_Gap2a_derep1_swarm.fasta",id=["Q.1","B.1.1.7"])

rule test3:
    input:
        expand(rez_dir + "/lineages/{id}/pair_id_Gap2a_derep1_swarm.fasta",id=["Q.1","B.1.1.7"]),
        #expand(rez_dir + "/lineages/{id}/pair_ref_Gap2a_derep1_swarm.fasta",id=["Q.1","B.1.1.7"])
        #rez_dir + "/lineages/Q.1/pair_ref_derep_data.csv"
        #rez_dir + "/pairwise_derep_results.csv" 

rule test4:
    input:
        ref = rez_dir + "/lineages/common/common_ref_derep_data.csv"
        #output_id = expand(rez_dir + "/lineages/{id}/common_id.fasta.gz", id = ["Q.1"])
        #expand(rez_dir + "/lineages/{id}/pair_id_Gap2a_derep1_swarm.fasta.swarminfo",id=["B.1.1.7"])


# expand(rez_dir + "/lineages/{id}/pair_id_Gap2a_derep1_swarm.fasta.swarminfo.cleaned_single.txt",id=["B.1.1.7"])
