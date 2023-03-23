import numpy as np
import os
import tqdm

cleanup_keep_fraction_of_columns = str(config["cleanup_keep_fraction_of_columns"])
cleanup_remove_proc_mostN_having = str(config["cleanup_remove_proc_mostN_having"])
reference_lineage = config["reference_lineage"]
lineages4summary = set(lineages) - set([config["reference_lineage"]])


rule do_subsample_per_periods:
    input:
        ids = rez_dir + "/lineages/{id}/common_id.txt",
        data = rez_dir + "/lineages/{id}/data.csv",
    output:
        partitions_info = rez_dir + "/lineages/{id}/periodic_partitions_F{perc,[0-9]+}S{seed,[0-9]+}_data.tsv",
        cluster = rez_dir + "/lineages/{id}/periodic_partitions/clF{perc,[0-9]+}S{seed,[0-9]+}"
    # notebook:
    #     "../notebooks/get_periodic_partitions.r.ipynb"
    script:
        "../notebooks/get_periodic_partitions.r.R"

rule do_subsample_per_periods_to_the_same_level:
    input:
        ids = rez_dir + "/lineages/{id}/common_id.txt",
        data = rez_dir + "/lineages/{id}/data.csv",
    output:
        partitions_info = rez_dir + "/lineages/{id}/periodic_partitions_N{numb,[0-9]+}S{seed,[0-9]+}_data.tsv",
        cluster = rez_dir + "/lineages/{id}/periodic_partitions/clN{numb,[0-9]+}S{seed,[0-9]+}"
    # notebook:
    #     "../notebooks/get_periodic_partitions_samelevel.r.ipynb"
    script:
        "../notebooks/get_periodic_partitions_samelevel.r.R"

rule subs1:
    input:
        expand(
            #rez_dir + "/lineages/{id}/periodic_partitions_F{perc}_S{seed}_data.tsv",
            rez_dir + "/lineages/{id}/periodic_partitions_analysis/clF{perc}S{seed}_veryfasttree.nw",
            perc = [2],
            seed = [2],
            id = lineages)


rule get_periodic_cluster_chronumental_tree:
    input:
        tree = rez_dir + "/lineages/{id}/periodic_partitions_analysis/clN{numb}S{seed}_veryfasttree_resolved.nw",
        metadata = rez_dir + "/lineages/{id}/data4timetree.csv"
    output:
        tree = rez_dir + "/lineages/{id}/periodic_partitions_analysis/clN{numb}S{seed}_chronumental.nw",
        node_data = rez_dir + "/lineages/{id}/periodic_partitions_analysis/clN{numb}S{seed}_chronumental_data.csv",
        calc_data = rez_dir + "/lineages/{id}/periodic_partitions_analysis/clN{numb}S{seed}_chronumental_data.txt",
    shell:
        """
        chronumental \
            --tree {input.tree} \
            --dates {input.metadata} \
            --tree_out {output.tree} \
            --dates_out {output.node_data} \
            --treat_mutation_units_as_normalised_to_genome_size 1 \
            --steps 10000 > {output.calc_data}
        """
rule parse_chronumentas_out_per_seed:
    input:
        calc_data = rez_dir + "/lineages/{id}/periodic_partitions_analysis/clN{numb}S{seed}_chronumental_data.txt",
    output:
        calc_data = rez_dir + "/lineages/{id}/periodic_partitions_analysis/clN{numb}S{seed}_chronumental_data_parsed.txt",
    conda:
        "../envs/R_env.yaml"
    notebook:
        "../notebooks/parse_chronumental_per_seed.r.ipynb"
        ""
rule subs2:
    input:
        expand(
            #rez_dir + "/lineages/{id}/periodic_partitions_F{perc}_S{seed}_data.tsv",
            #rez_dir + "/lineages/{id}/periodic_partitions_analysis/clN{numb}S{seed}_timetree_branch_lengths_er.txt",
            rez_dir + "/lineages/{id}/periodic_partitions_analysis/clN{numb}S{seed}_chronumental_data.txt",
            numb = [2000],
            seed = [2,1,3,4,5],
            id = lineages)


rule get_periodic_cluster_fasta:
    input:
        #fasta = rez_dir + "/lineages/{id}/alignment_nextclade.fasta.gz",
        fasta = rez_dir + "/lineages/{id}/common_id.fasta.gz",
        ids = rez_dir + "/lineages/{id}/periodic_partitions/cl{i}"
    output:
        rez_dir + "/lineages/{id}/periodic_partitions_analysis/cl{i, [^\_]+}.fasta"
    shell:
        """
            seqkit grep -w 0 -n -f {input.ids} {input.fasta} -o {output}
        """

rule get_periodic_cluster_tree1:
    input:
        rez_dir + "/lineages/{id}/periodic_partitions_analysis/cl{i}.fasta"
    output:
        rez_dir + "/lineages/{id}/periodic_partitions_analysis/cl{i, [^\_]+}_veryfasttree.nw"
    conda:
        "../envs/analysis.yaml"
    threads: 8
    shell:
        """
        cat {input} |  VeryFastTree -nosupport  -gamma -nt -gtr -out {output}  -double-precision  -threads {threads}
      """

rule get_periodic_cluster_resolved_tree:
    message:
        """
        Resolving multifurcatings
        """
    input:
        tree = rez_dir + "/lineages/{id}/periodic_partitions_analysis/cl{i}_veryfasttree.nw"
    output:
        tree = rez_dir + "/lineages/{id}/periodic_partitions_analysis/cl{i, [^\_]+}_veryfasttree_resolved.nw"
    conda:
        "../envs/analysis.yaml"
    threads: 1
    shell:
        """
        gotree resolve -i {input} -o {output}
        """

# rule reformatdata:
#     input:
#         metadata = rez_dir + "/lineages/{id}/data.csv"
#     output:
#         metadata = rez_dir + "/lineages/{id}/data4timetree.csv"
#     conda:
#         "../envs/R_env.yaml"
#     notebook:
#         "notebooks/preparedata4timetree.r.ipynb"
    

rule get_periodic_cluster_time_tree:
    """
    Refining tree
        - estimate timetree
        - use {params.coalescent} coalescent timescale
        - estimate {params.date_inference} node dates
        - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
    """
    input:
        tree = rez_dir + "/lineages/{id}/periodic_partitions_analysis/cl{i}_veryfasttree_resolved.nw",
        alignment = rez_dir + "/lineages/{id}/periodic_partitions_analysis/cl{i}.fasta",
        metadata = rez_dir + "/lineages/{id}/data4timetree.csv"
    output:
        tree = rez_dir + "/lineages/{id}/periodic_partitions_analysis/cl{i, [^\_]+}_timetree.nw",
        node_data = rez_dir + "/lineages/{id}/periodic_partitions_analysis/cl{i, [^\_]+}_timetree_branch_lengths.json"
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
rule extract_periodic_erate:
    input:
        rez_dir + "/lineages/{id}/periodic_partitions_analysis/cl{i, [^\_]+}_timetree_branch_lengths.json"
    output:
        rez_dir + "/lineages/{id}/periodic_partitions_analysis/cl{i, [^\_]+}_timetree_branch_lengths_er.txt"
    notebook:
        "../notebooks/extract_erate_from_json.py.ipynb" 

# rule analyse_periodic_subsamples:
#     input:
#         swarm = rez_dir + "/lineages/{id}/common_ids_peak_focus_Gap2a_derep1_swarm.fasta.swarminfo",
#         uc = rez_dir + "/lineages/{id}/common_ids_peak_focus_Gap2a_derep1.fasta.uc"
#     output:
#         rez = rez_dir + "/lineages/{id}/common_ids_peak_focus_Gap2a_derep1_swarm.data.txt"
#     conda:
#         "../envs/R_env.yaml"
#     notebook:
#         "../notebooks/peak_cluster_analysis.r.ipynb"


rule merge_data:
    input:
        expand(
            #rez_dir + "/lineages/{id}/periodic_partitions_F{perc}_S{seed}_data.tsv",
            rez_dir + "/lineages/{id}/periodic_partitions_analysis/clN{numb}S{seed}_timetree_branch_lengths_er.txt",
            numb = [2000,1000,500],
            seed = [2,1,3,4,5],
            id = lineages)
    output:
        rez_dir + "/erate_from_periodic_partition.txtx"
    shell:
        """
            cat {input} > {output}
        """

rule analyse_periodic_partitions_subsamples:
    input:
        rez_dir + "/erate_from_periodic_partition.txtx"
    output:
        rez_dir + "/erate_from_periodic_partition_erates.csv"
    log:
        notebook = rez_dir + "/erate_from_periodic_partition_erates.r.ipynb" 
    conda:
        "../envs/R_env.yaml"
    notebook:
        "../notebooks/erate_analysis4partitions.r.ipynb"