import numpy as np
import os
import tqdm

cleanup_keep_fraction_of_columns = str(config["cleanup_keep_fraction_of_columns"])
cleanup_remove_proc_mostN_having = str(config["cleanup_remove_proc_mostN_having"])
reference_lineage = config["reference_lineage"]
lineages4summary = set(lineages) - set([config["reference_lineage"]])


checkpoint do_subsample_per_periods_cluster:
    input:
        ids = rez_dir + "/lineages/{id}/common_id.txt",
        data = rez_dir + "/lineages/{id}/data.csv",
    output:
        cluster = directory(rez_dir + "/lineages/{id}/periodic_partitions_cluster/")
    notebook:
        "../notebooks/get_periodic_partitions_cluster.r.ipynb"

def aggregate_period_clusterr(wildcards):
    #checkpoint_output = checkpoints.choose_best_limit.get(**wildcards).output[1]
    checkpoint_output = checkpoints.do_subsample_per_periods_cluster.get(**wildcards).output[0]
    clsid = glob_wildcards(os.path.join(checkpoint_output, "week{i}.txt")).i
    out = expand(rez_dir + "/lineages/" + wildcards.id +"/periodic_partitions_cluster_analysis/week{i}_Gap2a_derep1_swarm.fasta",i=clsid)
    return(out)

def aggregate_uc(wildcards):
    #checkpoint_output = checkpoints.choose_best_limit.get(**wildcards).output[1]
    checkpoint_output = checkpoints.do_subsample_per_periods_cluster.get(**wildcards).output[0]
    clsid = glob_wildcards(os.path.join(checkpoint_output, "week{i}.txt")).i
    out = expand(rez_dir + "/lineages/" + wildcards.id +"/periodic_partitions_cluster_analysis/week{i}_Gap2a_derep1.fasta.uc",i=clsid)
    return(out)

def aggregate_swarm(wildcards):
    #checkpoint_output = checkpoints.choose_best_limit.get(**wildcards).output[1]
    checkpoint_output = checkpoints.do_subsample_per_periods_cluster.get(**wildcards).output[0]
    clsid = glob_wildcards(os.path.join(checkpoint_output, "week{i}.txt")).i
    out = expand(rez_dir + "/lineages/" + wildcards.id +"/periodic_partitions_cluster_analysis/week{i}_Gap2a_derep1_swarm.fasta.swarminfo",i=clsid)
    return(out)

rule get_periodic_cluster_fasta_4clustering:
    input:
        ids = rez_dir + "/lineages/{id}/periodic_partitions_cluster/week{i}.txt",
        fasta = rez_dir + "/lineages/{id}/common_id.fasta.gz",
        #fasta = rez_dir + "/lineages/{id}/alignment_nextclade.fasta.gz",
    output:
        rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/week{i, [^\_]+}.fasta.gz"
    shell:
        """
            seqkit grep -w 0 -n -f {input.ids} {input.fasta} -o {output}
        """


rule collect_clust_4_analysis:
    input:
        aggregate_period_clusterr
    output:
        rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all.fasta"
    shell:
        " cat {input} > {output} "

rule get_subset_originalALN:
    input:
        ids = rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_woSI.fasta.IDs",
        fasta = rez_dir + "/lineages/{id}/alignment_nextclade.fasta.gz",
    output:
        rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_matched_initial.fasta"
    shell:
        """
            seqkit grep -n -f {input.ids} {input.fasta} -o {output}
        """

rule get_periodic_cluster_tree1_clust:
    input:
        #rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_woSI.fasta"
        rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_matched_initial.fasta"
    output:
        rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_veryfasttree.nw"
    conda:
        "../envs/analysis.yaml"
    threads: 16
    shell:
        """
        cat {input} |  VeryFastTree -nosupport  -gamma -nt -gtr -out {output}  -double-precision  -threads {threads}
      """

rule get_periodic_cluster_resolved_tree_clust:
    message:
        """
        Resolving multifurcatings
        """
    input:
        tree = rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_veryfasttree.nw"
    output:
        tree = rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_veryfasttree_resolved.nw"
    conda:
        "../envs/analysis.yaml"
    threads: 1
    shell:
        """
        gotree resolve -i {input} -o {output}
        """
rule get_metadate4centroids:
    input:
        data = rez_dir + "/lineages/{id}/data.csv",
        ucs = aggregate_uc,
        swarm = aggregate_swarm, 
        cenroids = rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_woSI.fasta.IDs",
    output:
        metadata = rez_dir + "/lineages/{id}/data4timetree_4centroids.csv"
    conda:
        "../envs/R_env.yaml"
    notebook:
        "../notebooks/generate_metadata4centroid.r.ipynb"
        ""

rule get_periodic_cluster_time_tree_clust:
    """
    Refining tree
        - estimate timetree
        - use {params.coalescent} coalescent timescale
        - estimate {params.date_inference} node dates
        - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
    """
    input:
        tree = rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_veryfasttree_resolved.nw",
        alignment = rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_woSI.fasta",
        metadata = rez_dir + "/lineages/{id}/data4timetree.csv"
        #metadata = rez_dir + "/lineages/{id}/data4timetree_4centroids.csv"
    output:
        tree = rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_timetree.nw",
        node_data = rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_timetree_branch_length.json"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_filter_iqd = 2 # switched off now
        #clock_filter_iqd = 4 # switched off now
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

rule get_periodic_cluster_chronumental_tree_clust:
    """
    """
    input:
        #tree = rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_veryfasttree_resolved.nw",
        tree = rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_timetree.nw",
        alignment = rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_woSI.fasta",
        metadata = rez_dir + "/lineages/{id}/data4timetree.csv"
    output:
        tree = rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_chronumental.nw",
        node_data = rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_chronumental_data.csv", 
        calc_data = rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_chronumental_data.txt"
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
rule parse_chronumentas_out:
    input:
        calc_data = rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_chronumental_data.txt"
    output:
        calc_data = rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_chronumental_data_parsed.txt"
    conda:
        "../envs/R_env.yaml"
    notebook:
        "../notebooks/parse_chronumental.r.ipynb"
        ""


rule extract_erate_clust:
    input:
        rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_timetree_branch_length.json"
    output:
        rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_timetree_branch_length_er.txt"
    notebook:
        "../notebooks/extract_erate_from_json_cluswt.py.ipynb" 

rule aggregate_outputs:
    input:
        expand(rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_timetree_branch_length_er.txt", id = lineages)
    output:
        rez_dir + "/lineages/{id}/clusters4trees_analysis/merged_timetree_from_clust_branch_lengths_er.txt"
    shell:
        """
            cat {input} > {output}
        """

rule evaluate_tree_partition_clust:
    input:
        fasttree = expand(
            rez_dir + "/lineages/{id}/clusters4trees_analysis/merged_timetree_from_clust_branch_lengths_er.txt",
            #rez_dir + "/lineages/{id}/clusters_evals.csv",
            id = lineages  #["B.1.1.7"]
         ),
        chronumental = expand(
            rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_chronumental_data_parsed.txt",
            id = lineages
        )

    output:
            rez_dir + "/lineages/clusters4trees_clust_evolution_rate.csv",
    conda:
        "../envs/R_env.yaml"
    notebook:
        "notebooks/analyse_clusteredtreeER_clust.r.ipynb"

rule target_clust:
    input:
        rez_dir + "/lineages/clusters4trees_clust_evolution_rate.csv"

rule test_chron:
    input:
        expand(rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_chronumental.nw", id = lineages)
