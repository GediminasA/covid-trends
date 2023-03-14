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

rule get_periodic_cluster_fasta_4clustering:
    input:
        ids = rez_dir + "/lineages/{id}/periodic_partitions_cluster/week{i}.txt",
        fasta = rez_dir + "/lineages/{id}/common_id.fasta.gz",
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


rule get_periodic_cluster_tree1_clust:
    input:
        rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_woS.fasta"
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
        alignment = rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_woS.fasta",
        metadata = rez_dir + "/lineages/{id}/data4timetree.csv"
    output:
        tree = rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_timetree.nw",
        node_data = rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_timetree_branch_length.json"
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
rule targettt_clust:
    input:
        expand(rez_dir + "/lineages/{id}/periodic_partitions_cluster_analysis/all_timetree.nw", id = lineages)

