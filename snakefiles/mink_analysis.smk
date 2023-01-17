outdir = config["work_dir"]+"/mink_phylogeny"
#lineage = lineage.split()
rule minks:
    input:
        #data = rez_dir + "/lineages/mink/data.csv",
        rez = rez_dir + "/initial_set_minimap2_match_all.sam"
        #fasta = rez_dir+"/tmp/lineages_all.fasta"

rule get_all_clean_ids:
    input:
        nextclade = rez_dir+"/nextclade_report.csv"
    output:
        fasta = rez_dir+"/tmp/lineages_all.id"
    threads: 12
    notebook:
        "notebooks/filter_out_bad.r.ipynb"



rule get_all_clean_fasta:
    input:
        ids = rez_dir+"/tmp/lineages_all.id"
    output:
        fasta = rez_dir+"/tmp/lineages_all.fasta"
    params:
        alignment = config["sequences"],
    threads: 12
    shell:
        """
            seqkit grep -w 0 -n -f {input.ids} {params.alignment} -o {output.fasta}
        """
rule get_animal_ids:
    input:
        pangolin = rez_dir+"/pangolin_lineage_report.csv",
        nextclade = rez_dir+"/nextclade_report.csv",
        meta = config["meta"]
    output:
        ids = rez_dir + "/lineages/mink/ids.txt",
        data = rez_dir + "/lineages/mink/data.csv",
    params:
        animal = "mink"
    threads: 12
    notebook:
        "notebooks/filter_out_animal.r.ipynb"

rule extract_animal_seqs:
    input:
        ids = rez_dir + "/lineages/mink/ids.txt",
    output:
        aln = rez_dir + "/lineages/mink/sequences.fasta",
    params:
        alignment = config["sequences"],
    threads: 12
    shell:
        """
            seqkit grep -w 0 -n -f {input.ids} {params.alignment} -o {output.aln}
        """
   

rule get_hisat_index:
    input: "{stem}.fasta"
    output: "{stem}_hisat2.1.ht2"
    params: stem = "{stem}_hisat2"
    threads: 80
    conda:
        "../envs/searchsimilar.yaml"
    shell: "hisat2-build -p {threads} {input} {params.stem}"    

rule remove_n_subs:
    input:
        seqs = rez_dir + "/lineages/mink/sequences.fasta",
    output:
        seqs = rez_dir + "/lineages/mink/sequences_woN.fasta",
    conda:
        "../envs/searchsimilar.yaml"
    shell:
       " seqkit seq -g -G N -i {input} -o {output} "

rule remove_n_alls:
    input:
        fasta = config["sequences"]
    output:
        raw_woN = temp(rez_dir+"/tmp/lineages_all_woN.fasta")
    conda:
        "../envs/searchsimilar.yaml"
    shell:
       " seqkit seq -g -G N -i {input} -o {output} "

rule closestv2:
    input:
        reference = rez_dir + "/lineages/mink/sequences_woN.fasta",
        reference_index = rez_dir + "/lineages/mink/sequences_woN_hisat2.1.ht2",
        fasta = rez_dir+"/tmp/lineages_all_woN.fasta"
    output:
        rez = rez_dir+"/initial_set_minimap2_match_all.sam"
    params:
        index = rez_dir + "/lineages/mink/sequences_woN_hisat2",
    conda:
        "../envs/searchsimilar.yaml"
    threads: 80
    shell:
        """
        hisat2 -p {threads}  -x {params.index}  -f {input.fasta}   --no-spliced-alignment  -k 2 > {output}
        """

rule sam_paf:
    input: "{stem}.sam"
    output: "{stem}.paf"
    conda:
        "envs/analysis.yaml"
    shell: " paftools.js sam2paf -p  {input} > {output}"

rule sortminimap:
    input:
        rez2 = rez_dir+"/initial_set_minimap2_match_all.paf"
    output: 
        ids = outdir+"/tmp/minimap2_match_chosen.ids",
        ids_target = outdir+"/tmp/minimap2_match_chosen.ids_target",
        ids_additional = outdir+"/tmp/minimap2_match_chosen.ids_additional",
    conda:
        "envs/r.yaml"
    params:
        target_amount = config["target_amount"] 
    script:
        "scripts/minimap.R"    


# below borrowed rules for pphylogenetic analyses
rule extract_top_seqyences:
    input:
        ids = outdir+"/tmp/minimap2_match_chosen.ids",
        fasta = rez_dir+"/tmp/lineages_all.fasta"
    output:
       fasta_ini = outdir+"/tmp/minimap2_match_chosen_ini_not_ordered.fa",
    threads: 80
    conda:
        "envs/analysis.yaml"
    shell:
        "seqkit grep -n -j {threads} -f {input.ids} -o {output} -w 0 -i {input.fasta} "


rule get_top_sequences:
    input:
       ids = outdir+"/tmp/minimap2_match_chosen.ids",
       fasta_ini = outdir+"/tmp/minimap2_match_chosen_ini_not_ordered.fa",
    output:
       fasta_final = outdir+"/tmp/minimap2_match_chosen.fasta",
    conda:
        "envs/analysis.yaml"
    shell:
       "  rm -f {input.fasta_ini}.fai ; seqkit faidx {input.fasta_ini} --infile-list {input.ids} -w 0 -o {output} "



rule add_sub_sequences:
    input:
        fasta_chosen = outdir+"/tmp/minimap2_match_chosen.fasta",
        fasta_sub = outdir+"/tmp/lineages_sub.fasta",
    output:
        sequences = outdir+"/tmp/minimap2_match_chosen_plus_sub.fasta",
    shell:
        " cat {input} > {output} "


rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = outdir+"/tmp/minimap2_match_chosen_plus_sub.fasta",
        reference = config["reference"]
    output:
        alignment = outdir+"/tmp/minimap2_match_chosen_aligned.fasta",
        #alignment = outdir+"/tmp/minimap2_match_chosen_{lineage,\D{1,3}?\.\d+\.\d+\.\d+}_aligned.fasta",
    conda:
        "envs/analysis.yaml"
    threads: 20
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --fill-gaps \
            --remove-reference \
            --nthreads {threads}
        """


rule inital_tree:
    message:
        """
        Running very fast tree
        """
    input:
        alignment = outdir+"/tmp/minimap2_match_chosen_aligned.fasta",
    output:
        alignment = outdir+"/tmp/minimap2_match_chosen_aligned_veryfast.nwk",
    conda:
        "envs/analysis.yaml"
    threads: 20
    shell:
        """
        cat {input} |  VeryFastTree -nosupport  -gamma -nt -gtr -out {output}  -double-precision  -threads {threads}
      """
        # export OMP_NUM_THREADS={threads} ; cat {input} |  FastTree -nosupport  -gamma -nt -gtr -out {output}  
        #cat {input} |  VeryFastTree -nosupport  -gamma -nt -gtr -out {output}  -threads {threads}



rule resolve_tree:
    message:
        """
        Resilving multifurcatings
        """
    input:
        tree = outdir+"/tmp/minimap2_match_chosen_aligned_veryfast.nwk",
    output:
        tree = outdir+"/tmp/minimap2_match_chosen_aligned_veryfast_resolved.nwk",
    conda:
        "envs/analysis.yaml"
    threads: 1
    shell:
        """
        gotree resolve -i {input} -o {output}
        """

rule fix_brancges:
    message:
        """
        Recalculate branch length
        """
    input:
        tree = outdir+"/tmp/minimap2_match_chosen_aligned_veryfast_resolved.nwk",
        alignment = outdir+"/tmp/minimap2_match_chosen_aligned.fasta",
    output:
        tree = outdir+"/output/initialtree.nwk",
    conda:
        "envs/raxml.yaml"
    params:
        prefix = outdir+"/tmp/minimap2_match_chosen_aligned_veryfast_resolved_fixed_raxml",
        output = outdir+"/tmp/minimap2_match_chosen_aligned_veryfast_resolved_fixed_raxml.raxml.bestTree"
    threads: 20
    shell:
        """
        raxml-ng --blopt nr_safe  --redo --precision 12--redo --msa {input.alignment} --model GTR+FO+I+R4 --evaluate --tree {input.tree} --prefix {params.prefix} --threads auto
        mv {params.output} {output}
        """

rule get_time_tree:
    """
    Refining tree
        - estimate timetree
        - use {params.coalescent} coalescent timescale
        - estimate {params.date_inference} node dates
        - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
    """
    input:
        tree = outdir+"/output/initialtree.nwk",
        alignment = outdir+"/tmp/minimap2_match_chosen_aligned.fasta",
        metadata = config["meta"],
    output:
        tree = outdir+"/output/timetree.nwk",
        node_data = outdir + "/output/timetree_branch_lengths.json"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_filter_iqd = 4
    conda:
        "envs/analysis.yaml"
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
            --clock-filter-iqd {params.clock_filter_iqd} \
        """


rule get_outgroup:
    input:
        timetree = outdir+"/output/timetree.nwk",
    output:
        timetreeroot = outdir+"/tmp/timetree.root",
    shell:
        """
        ./scripts/extract_outgroup_tree.jl -i {input} -o {output}
        """



rule get_matching_original_tree:
    input:
        tree = outdir+"/output/initialtree.nwk",
        timetree = outdir+"/output/timetree.nwk",
    output:
        treeajusted = outdir+"/tmp/initialtreeTrimmed.nwk",
    conda:
        "envs/analysis.yaml"
    shell:
        """
        gotree prune -i {input.tree} -c {input.timetree} -o {output.treeajusted}
        """


rule get_reroted_matching_original_tree:
    input:
        treeajusted = outdir+"/tmp/initialtreeTrimmed.nwk",
        timetreeroot = outdir+"/tmp/timetree.root",
    output:
        treeajusted = outdir+"/tmp/initialtreeTrimmedRooted.nwk",
    conda:
        "envs/analysis.yaml"
    shell:
        """
        gotree reroot outgroup -l {input.timetreeroot} -i {input.treeajusted} -o {output}
        """

rule split_tree:
    input:
        treeajusted = outdir+"/tmp/initialtreeTrimmedRooted.nwk",
    output:
        directory(outdir+"/tmp/initialtreeTrimmedRooted_splits"),
    params:
        maxn = config["maximum_amount_4_splits"],
    conda:
        "envs/analysis.yaml"
    shell:
        """
        ./scripts/split_tree.jl -i {input} -o {output}
        """



rule get_clustering:
    input:
        outdir+"/tmp/initialtreeTrimmedRooted_splits",
    output:
        outdir+"/output/clusters.txt",
    params:
        wd = outdir+"/tmp/initialtreeTrimmedRooted_splits"
    conda:
        "envs/phydelity.yaml"
    shell:
        """
        wd=` pwd `
        cd {params.wd}
        for l in tree*.nwk ; do
            stem=${{l/\.nwk/""}} 
            echo $stem
            nohup phydelity.py --k 2 --collapse_zero_branch_length -t  $l > pydelity_wd_$l.log &
        done
        wait
        for l in tree*.nwk ; do
            stem=${{l/\.nwk/""}} 
            clsf=cluster_phydelity_k2_sol0_"$stem".txt
            cat $clsf | grep -v "CLUSTER" | sed "s/^/cl"$stem"_/g" >> $wd"/"{output}
        done
        """

rule get_focused_with_siblings:
    input:
        treeajusted = outdir+"/output/timetree.nwk",
        lineage_seqs_sub_ids = outdir+"/tmp/lineage_sub.ids.txt",
    output:
        treefocused = outdir+"/output/focused.nwk",
    conda:
        "envs/analysis.yaml"
    shell:
        """
        ./scripts/focus_tree.jl -i {input.treeajusted} -l {input.lineage_seqs_sub_ids} -o {output} -s
        """

rule get_focused_wo_siblings:
    input:
        treeajusted = outdir+"/output/timetree.nwk",
        lineage_seqs_sub_ids = outdir+"/tmp/lineage_sub.ids.txt",
    output:
        treefocused = outdir+"/output/focused_wosiblings.nwk",
    conda:
        "envs/analysis.yaml"
    shell:
        """
        ./scripts/focus_tree.jl -i {input.treeajusted} -l {input.lineage_seqs_sub_ids} -o {output} 
        """

rule traits:
    message:
        """
        Inferring ancestral traits for {params.columns!s}
          - increase uncertainty of reconstruction by {params.sampling_bias_correction} to partially account for sampling bias
        """
    input:
        tree = outdir+"/output/timetree.nwk",
        metadata = config["meta"],
    output:
        node_data = outdir + "/output/timetree_traits.json"
    params:
        columns = "country",
        #columns = "region country",
        sampling_bias_correction = 3
    conda:
        "envs/analysis.yaml"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence \
            --sampling-bias-correction {params.sampling_bias_correction}
        """

rule collapse_tree:
    input:
        node_data = outdir + "/output/timetree_traits.json",
        treefocused = outdir+"/output/focused_wosiblings.nwk",
    output:
        treefocusedcollapsed = outdir+"/output/focused_wosiblings_collapsed1.nwk",
        treefocusedcollapsed_data = outdir+"/output/focused_wosiblings_collapsed1.json",
    conda:
        "envs/analysis.yaml"
    shell:
        """
        ./scripts/collapse_tree.jl -j {input.node_data} -i {input.treefocused} -o {output.treefocusedcollapsed} -d {output.treefocusedcollapsed_data} 
        """

rule collapse_tree2:
    input:
        node_data = outdir + "/output/timetree_traits.json",
        treefocused = outdir+"/output/focused.nwk",
    output:
        treefocusedcollapsed = outdir+"/output/focused_collapsed2.nwk",
        treefocusedcollapsed_data = outdir+"/output/focused_collapsed2.json", 
    conda:
        "envs/analysis.yaml"
    shell:
        """
        ./scripts/collapse_tree.jl -j {input.node_data} -i {input.treefocused} -o {output.treefocusedcollapsed} -d {output.treefocusedcollapsed_data} 
        """


# rule extract_top_seqyences:
#     input:
#         ids = outdir+"/tmp/minimap2_match_chosen.ids",
#         fasta = outdir+"/tmp/lineages_all.fasta"
#     output:
#        fasta_ini = outdir+"/tmp/minimap2_match_chosen_ini_not_ordered.fa",
#     threads: 80
#     conda:
#         "envs/analysis.yaml"
#     shell:
#         "seqkit grep -n -j {threads} -f {input.ids} -o {output} -w 0 -i {input.fasta} "


# rule get_top_sequences:
#     input:
#        ids = outdir+"/tmp/minimap2_match_chosen.ids",
#        fasta_ini = outdir+"/tmp/minimap2_match_chosen_ini_not_ordered.fa",
#     output:
#        fasta_final = outdir+"/tmp/minimap2_match_chosen.fasta",
#     conda:
#         "envs/analysis.yaml"
#     shell:
#        "  rm -f {input.fasta_ini}.fai ; seqkit faidx {input.fasta_ini} --infile-list {input.ids} -w 0 -o {output} "


       