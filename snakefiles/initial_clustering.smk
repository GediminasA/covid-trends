import numpy as np
import os
import tqdm

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
        rows_after_filter = rez_dir + "/lineages/{id}/alignment_nextclade_filtered_keptrows_{k}_{t}.txt", 
        columns_after_filter = rez_dir + "/lineages/{id}/alignment_nextclade_filtered_keptcolumns_{k}_{t}.txt", 
    params:
        keepPositionsFraction = "{k}",#0.8, # 0.7 0.8 0.9 0.95
        topNhaving_sequences_to_remove ="{t}", # 10, # 10 20 10 5 2
        reference_seq_id = "MN908947",
        country = "Lithuania",
    threads: 1000
    log: rez_dir + "/lineages/{id}/alignment_nextclade_filtered_parsed_{k}_{t}.ipynb"
    shell:
        """
        papermill  scripts/julia_modules/JuliaClusterAndTreeTools/notebooks/filter_aln.ipynb \
        {log} \
        -p in_file_name {input.in_file_name} \
        -p rows_after_filter {output.rows_after_filter} \
        -p columns_after_filter {output.columns_after_filter} \
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

rule remove_size_info:
    input:
        "{stem}.fasta"
    output:
        "{stem}_woS.fasta"
    threads: 12
    conda: "../envs/clustering_tools.yaml"
    shell:
        """
            vsearch --fastx_filter {input} --fastaout {output} --xsize
            sed -i "s/'//g" {output}
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


rule swarmD2:
    input:
        "{stem}.fasta",
    output:
        fa = "{stem}_swarmD2.fasta",
        uc = "{stem}_swarmD2.fasta.swarminfo",
        stru = "{stem}_swarmD2.fasta.internstr"
    params:
    conda: "../envs/clustering_tools.yaml"
    threads: 1000
    shell:
        '''
        swarm  --threads {threads} -z  -d  2   -w {output[0]} -r -i {output.stru} -o {output.uc} {input}
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

rule swarmNF: #NF mean non fastidious
    input:
        "{stem}.fasta",
    output:
        fa = "{stem}_swarmNF.fasta", 
        uc = "{stem}_swarmNF.fasta.swarminfo",
        stru = "{stem}_swarmNF.fasta.internstr"
    params:
    conda: "../envs/clustering_tools.yaml"
    threads: 1000
    shell:
        '''
        swarm --threads {threads} -z -d 1   -w {output[0]} -r -i {output.stru} -o {output.uc} {input}
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
    threads: 1000
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
        ids = expand(rez_dir + "/lineages/{id}/common_id_Gap2a_derep1.fasta.uc", id = lineages),
        setupmark = config["work_dir"]+"/RsetupDone.txt"
    log:
        notebook = rez_dir + "/lineages/common/common_ref_derep_data.r.ipynb"
    output:
        ref = rez_dir + "/lineages/common/common_ref_derep_data.csv"
    conda:
        "../envs/R_env.yaml"
    notebook:
        "../notebooks/analyse_derep_common.r.ipynb"

rule amalyze_common_derep_swarm:
    input:
        datas = expand(rez_dir + "/lineages/{id}/data.csv", id = lineages),
        ids = expand(rez_dir + "/lineages/{id}/common_id_Gap2a_derep1.fasta.uc", id = lineages),
        ids2 = expand(rez_dir + "/lineages/{id}/common_id_Gap2a_derep1_swarm.fasta.internstr", id = lineages),
        setupmark = config["work_dir"]+"/RsetupDone.txt"
    log:
        notebook = rez_dir + "/lineages/common/common_ref_derep_swarm_data.r.ipynb"
    output:
        ref = rez_dir + "/lineages/common/common_ref_derep_swarm_data.csv"
    conda:
        "../envs/R_env.yaml"
    notebook:
        "../notebooks/analyse_derep_swarm_common.r.ipynb"

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

        
# test simple evolution rate inference

rule get_peak_ids:
    input:
        ids = rez_dir + "/lineages/{id}/common_id.txt",
        data = rez_dir + "/lineages/{id}/data.csv",
    output:
        ids_focus = rez_dir + "/lineages/{id}/common_ids_peak_focus.txt",
        ids_other = rez_dir + "/lineages/{id}/common_ids_peak_other.txt",
    notebook:
        "../notebooks/get_peak_ids.r.ipynb"

rule analyse_peaks:
    input:
        swarm = rez_dir + "/lineages/{id}/common_ids_peak_focus_Gap2a_derep1_swarm.fasta.swarminfo",
        uc = rez_dir + "/lineages/{id}/common_ids_peak_focus_Gap2a_derep1.fasta.uc"
    output:
        rez = rez_dir + "/lineages/{id}/common_ids_peak_focus_Gap2a_derep1_swarm.data.txt"
    conda:
        "../envs/R_env.yaml"
    notebook:
        "../notebooks/peak_cluster_analysis.r.ipynb"


rule get_peak_ids_fasta_focus:
    input:
        ids = rez_dir + "/lineages/{id}/common_ids_peak_focus.txt",
        fasta = rez_dir + "/lineages/{id}/common_id.fasta.gz"
    output:
        fasta = rez_dir + "/lineages/{id}/common_ids_peak_focus.fasta.gz",
    shell:
        """
            seqkit grep -n -f {input.ids} {input.fasta} -o {output.fasta}
        """

rule print_common_ids:
    input:
        fasta = rez_dir + "/lineages/{id}/common_id.fasta.gz"
    output:
        ids = rez_dir + "/lineages/{id}/common_id.txt"
    shell:
        """
            seqkit seq  -n -i  {input} -w 0 -o {output}
        """

rule decenttree:
    input:
        #rez_dir + "/lineages/{id}/common_id_Gap2a_derep1.fasta"
        rez_dir + "/lineages/{id}/common_id_Gap2a_derep1_swarmNF_woS.fasta"
    output:
        rez_dir + "/lineages/{id}/common_id_decenttree.nwk"
    threads: 36
    shell:
        """
            external_programs/decenttree/bin/decenttree -nt {threads} \
                -fasta {input} -max-dist 1000	-t NJ-R-V \
                -out {output} -bar
        """


rule clustertree:
    input:
        tree = rez_dir + "/lineages/{id}/common_id_decenttree.nwk"
    output:
        cluster = rez_dir + "/lineages/{id}/clusters/{p}_ini.tsv",
    shell:
        "TreeCluster.py  -i {input.tree} -t {wildcards.p} -m length_clade  -o {output} "
        #"TreeCluster.py  -i {input.tree} -t {wildcards.p} -m single_linkage_union -o {output} "

#single_linkage 1
#single_linkage_cut 1
#single_linkage_union 1
#sum_branch 0
#sum_branch_clade 0
#med_clade 1
#root_dist 1
#avg_clade 1
#leaf_dist_max 1
#length 1
#length_clade 2
#max 1

#avg_clade, leaf_dist_avg, leaf_dist_max, leaf_dist_min, length, length_clade, max, max_clade,
#med_clade, root_dist, single_linkage, single_linkage_cut, single_linkage_union, sum_branch, sum_branch_clade

# rule clustertrees:
#     input:
#         [rez_dir + "/lineages/{id}/clusters/"+str(p)+".tsv" for p in np.linspace(0.0001, 0.0005, 1000)]
#     output:
#         cluster = rez_dir + "/lineages/{id}/clusters/generated.txt"
#     shell:
#         "touch {output}"


rule merge_cluster4evaluatecluster:
    input:
        swarm = rez_dir + "/lineages/{id}/common_id_Gap2a_derep1_swarmNF.fasta.swarminfo",
        darep = rez_dir + "/lineages/{id}/common_id_Gap2a_derep1.fasta.uc",
    output:
        merged = rez_dir + "/lineages/{id}/common_id_Gap2a_derep1_swarmNF_mergedclusters.csv",
    conda:
        "../envs/R_env.yaml"
    notebook:
        "../notebooks/join_clusters.r.ipynb"

rule evaluatecluster:
    input:
        cluster = rez_dir + "/lineages/{id}/clusters/{p}_ini.tsv",
        data = rez_dir + "/lineages/{id}/data.csv",
        merged = rez_dir + "/lineages/{id}/common_id_Gap2a_derep1_swarmNF_mergedclusters.csv",
    output:
        rez_dir + "/lineages/{id}/clusters/{p}_eval.tsv"
    conda:
        "../envs/R_env.yaml"
    script:
       "../notebooks/evaluate_cluster.r.R"
    # notebook:
    #     "../notebooks/evaluate_cluster.r.ipynb"

rule merge_evaluations:
    input:
        [rez_dir + "/lineages/{id}/clusters/"+str(p)+"_eval.tsv" for p in np.arange(0.00001, 0.0005, 0.00001)]
    output:
        rez_dir + "/lineages/{id}/clusters_evals.csv"
    shell:
        "cat {input} > {output}"


checkpoint choose_best_limit:
    input:
        evals = rez_dir + "/lineages/{id}/clusters_evals.csv",
        data = rez_dir + "/lineages/{id}/data.csv",
        merged = rez_dir + "/lineages/{id}/common_id_Gap2a_derep1_swarmNF_mergedclusters.csv",
    output:
        clusters = directory(rez_dir + "/lineages/{id}/clusters4trees"), 
    params:
        clusters = rez_dir + "/lineages/{id}/clusters/",
        max_number_of_cluster_to_write = 10
    log:
        notebook = rez_dir + "/lineages/{id}/clusters_evals_best.ipynb",
    conda:
        "../envs/R_env.yaml"
    notebook:
        "../notebooks/clusters_evals_best_choose.r.ipynb"

def aggregate_input(wildcards):
    #checkpoint_output = checkpoints.choose_best_limit.get(**wildcards).output[1]
    checkpoint_output = checkpoints.choose_best_limit.get(**wildcards).output[0]
    clsid = glob_wildcards(os.path.join(checkpoint_output, "cl{i}")).i
    #out = expand(rez_dir + "/lineages/"+wildcards.id+"/clusters4trees_analysis/cl{i}.fasta",i=clsid)
    #out = expand(rez_dir + "/lineages/"+wildcards.id+"/clusters4trees_analysis/cl{i}_veryfasttree_resolved.nw",i=clsid)
    #out = expand(rez_dir + "/lineages/"+wildcards.id+"/clusters4trees_analysis/cl{i}_timetree.nw",i=clsid)
    out = expand(rez_dir + "/lineages/"+wildcards.id+"/clusters4trees_analysis/cl{i}_timetree_branch_lengths_er.txt",i=clsid)
    #import sys 
    #sys.exit()
    return(out)

rule get_cluster_fasta:
    input:
        fasta = rez_dir + "/lineages/{id}/alignment_nextclade.fasta.gz",
        ids = rez_dir + "/lineages/{id}/clusters4trees/cl{i}"
    output:
        rez_dir + "/lineages/{id}/clusters4trees_analysis/cl{i}.fasta"
    shell:
        """
            seqkit grep -w 0 -n -f {input.ids} {input.fasta} -o {output}
        """

rule get_cluster_tree1:
    input:
        rez_dir + "/lineages/{id}/clusters4trees_analysis/cl{i}.fasta"
    output:
        rez_dir + "/lineages/{id}/clusters4trees_analysis/cl{i}_veryfasttree.nw"
    conda:
        "../envs/analysis.yaml"
    threads: 32
    shell:
        """
        cat {input} |  VeryFastTree -nosupport  -gamma -nt -gtr -out {output}  -double-precision  -threads {threads}
      """

rule get_cluster_resolved_tree:
    message:
        """
        Resolving multifurcatings
        """
    input:
        tree = rez_dir + "/lineages/{id}/clusters4trees_analysis/cl{i}_veryfasttree.nw"
    output:
        tree = rez_dir + "/lineages/{id}/clusters4trees_analysis/cl{i}_veryfasttree_resolved.nw"
    conda:
        "../envs/analysis.yaml"
    threads: 1
    shell:
        """
        gotree resolve -i {input} -o {output}
        """

rule reformatdata:
    input:
        metadata = rez_dir + "/lineages/{id}/data.csv"
    output:
        metadata = rez_dir + "/lineages/{id}/data4timetree.csv"
    conda:
        "../envs/R_env.yaml"
    notebook:
        "notebooks/preparedata4timetree.r.ipynb"
    

rule get_cluster_time_tree:
    """
    Refining tree
        - estimate timetree
        - use {params.coalescent} coalescent timescale
        - estimate {params.date_inference} node dates
        - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
    """
    input:
        tree = rez_dir + "/lineages/{id}/clusters4trees_analysis/cl{i}_veryfasttree_resolved.nw",
        alignment = rez_dir + "/lineages/{id}/clusters4trees_analysis/cl{i}.fasta",
        metadata = rez_dir + "/lineages/{id}/data4timetree.csv"
    output:
        tree = rez_dir + "/lineages/{id}/clusters4trees_analysis/cl{i}_timetree.nw",
        node_data = rez_dir + "/lineages/{id}/clusters4trees_analysis/cl{i}_timetree_branch_lengths.json"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_filter_iqd = 4 # switched off now
    conda:
        "../envs/analysis.yaml"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            --date-confidence  --root best \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

rule extract_erate:
    input:
        rez_dir + "/lineages/{id}/clusters4trees_analysis/cl{i}_timetree_branch_lengths.json"
    output:
        rez_dir + "/lineages/{id}/clusters4trees_analysis/cl{i}_timetree_branch_lengths_er.txt"
    notebook:
        "../notebooks/extract_erate_from_json.py.ipynb" 


rule aggregate:
    input:
        aggregate_input
    output:
        rez_dir + "/lineages/{id}/clusters4trees_analysis/merged_timetree_branch_lengths_er.txt"
    shell:
        """
            touch {output}
            cat {input} > {output}
        """


rule evaluate_tree_partition:
    input:
        expand(
            rez_dir + "/lineages/{id}/clusters4trees_analysis/merged_timetree_branch_lengths_er.txt",
            #rez_dir + "/lineages/{id}/clusters_evals.csv",
            id = lineages  #["B.1.1.7"]
         )
    output:
            rez_dir + "/lineages/clusters4trees_evolution_rate.csv",
    conda:
        "../envs/R_env.yaml"
    notebook:
        "notebooks/analyse_clusteredtreeER.r.ipynb"



rule test2:
    input:
        #expand(rez_dir + "/lineages/{id}/kept_rows.txt",id=["BA.2"])
        expand(rez_dir + "/lineages/{id}/pair_id_Gap2a_derep1_swarm.fasta",id=["Q.1","B.1.1.7"]),
        expand(rez_dir + "/lineages/{id}/pair_ref_Gap2a_derep1_swarm.fasta",id=["Q.1","B.1.1.7"])

rule test3:
    input:
        #expand(rez_dir + "/lineages/{id}/pair_id_Gap2a_derep1_swarm.fasta",id=["Q.1","B.1.1.7"]),
        #expand(rez_dir + "/lineages/{id}/pair_ref_Gap2a_derep1_swarm.fasta",id=["Q.1","B.1.1.7"])
        #rez_dir + "/lineages/Q.1/pair_ref_derep_data.csv"
        rez_dir + "/pairwise_derep_results.csv" 

rule test4:
    input:
        ref = rez_dir + "/lineages/common/common_ref_derep_data.csv"
        #output_id = expand(rez_dir + "/lineages/{id}/common_id.fasta.gz", id = ["Q.1"])
        #expand(rez_dir + "/lineages/{id}/pair_id_Gap2a_derep1_swarm.fasta.swarminfo",id=["B.1.1.7"])


rule test5:
    input:
        #ref = rez_dir + "/lineages/common/common_ref_derep_data.csv"
        #output_id = expand(rez_dir + "/lineages/{id}/common_id.fasta.gz", id = ["Q.1"])
        #expand(rez_dir + "/lineages/{id}/common_id_Gap2a_derep1_swarm.fasta.internstr", id = lineages),
        #ref = rez_dir + "/lineages/common/common_ref_derep_swarm_data.csv"
        #expand(rez_dir + "/lineages/{id}/common_ids_peak.txt",id=lineages)
        # "audines11/rez/lineages/BA.2.9/common_ids_peak_focus_Gap2a_derep1_swarm.data.txt"
        expand(rez_dir + "/lineages/{id}/common_ids_peak_focus_Gap2a_derep1_swarm.data.txt", id = lineages)

rule test6:
    input:
        rez_dir + "/lineages/clusters4trees_evolution_rate.csv",
# expand(rez_dir + "/lineages/{id}/pair_id_Gap2a_derep1_swarm.fasta.swarminfo.cleaned_single.txt",id=["B.1.1.7"])
