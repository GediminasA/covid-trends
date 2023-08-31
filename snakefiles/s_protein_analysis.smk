cleanup_keep_fraction_of_columns = str(config["cleanup_keep_fraction_of_columns"])
cleanup_remove_proc_mostN_having = str(config["cleanup_remove_proc_mostN_having"])
reference_lineage = config["reference_lineage"]
lineages4summary = set(lineages) - set([config["reference_lineage"]])

rule extract_reference_sequence_s_protein:
    input:
        alignment = rez_dir+"/gene_S.translation.fasta.gz" 
    output:
        ref = rez_dir + "/lineages/alignment_nextclade_refseq_sprot.fasta.gz",
    params:
        id = config["reference_id"]
    threads: 12
    conda: 
        "../envs/clustering_tools.yaml"
    shell:
        "seqkit grep -w 0 -n -p '{params.id}' {input.alignment} -o {output.ref}  "

rule extract_alignment_s_protein:
    input:
        ids = rez_dir + "/lineages/{id}/ids.txt",
        alignment = rez_dir+"/gene_S.translation.fasta.gz", 
        ref = rez_dir + "/lineages/alignment_nextclade_refseq_sprot.fasta.gz",
    output:
        aln = rez_dir + "/lineages/{id}/sprotein/all_sequences.fasta.gz",
    threads: 12
    conda: 
        "../envs/clustering_tools.yaml"
    shell:
        """  
            seqkit grep -w 0 -n -f {input.ids} {input.alignment} -o {output.aln}
            cat {input.ref} >> {output.aln}
        """
rule extract_alignment_s_protein_lt:
    input:
        ids = rez_dir + "/lineages/{id}/ids_lt.txt",
        alignment = rez_dir+"/gene_S.translation.fasta.gz", 
        ref = rez_dir + "/lineages/alignment_nextclade_refseq_sprot.fasta.gz",
    output:
        aln = rez_dir + "/lineages/{id}/sprotein/lt_sequences.fasta.gz",
    threads: 12
    conda: 
        "../envs/clustering_tools.yaml"
    shell:
        """  
            seqkit grep -w 0 -n -f {input.ids} {input.alignment} -o {output.aln}
            cat {input.ref} >> {output.aln}
        """
rule analyse_conservation:
    input:
        aln = rez_dir + "/lineages/{id}/sprotein/all_sequences.fasta.gz"
    output:
        csv = rez_dir + "/lineages/{id}/sprotein/all_{id}_sequences_conservation_all.csv"
    log:
        log = rez_dir + "/lineages/{id}/sprotein/all_{id}_sequences_conservation_all.ipynb"
    threads: 16
    conda:
        "../envs/R_env.yaml"
    shell:
        """
        papermill scripts/julia_modules/JuliaClusterAndTreeTools/notebooks/sprot_entropy.ipynb \
        {log} \
        -p fasta {input.aln} \
        -p  output {output.csv} 
        """

rule analyse_conservation_lt:
    input:
        aln = rez_dir + "/lineages/{id}/sprotein/lt_sequences.fasta.gz"
    output:
        csv = rez_dir + "/lineages/{id}/sprotein/lt_{id}_sequences_conservation_all.csv"
    log:
        log = rez_dir + "/lineages/{id}/sprotein/lt_{id}_sequences_conservation_all.ipynb"
    threads: 16
    conda:
        "../envs/R_env.yaml"
    shell:
        """
        papermill scripts/julia_modules/JuliaClusterAndTreeTools/notebooks/sprot_entropy.ipynb \
        {log} \
        -p fasta {input.aln} \
        -p  output {output.csv} 
        """

rule analyse_reference_based:
    input:
        data = rez_dir + "/lineages/{id}/data.csv",
    output:
        csv = rez_dir + "/lineages/{id}/sprotein/reference_based_mutations_{id}_conservation.csv",
        data = rez_dir + "/lineages/{id}/sprotein/reference_based_mutations_{id}_data.csv"
    log:
        log = rez_dir + "/lineages/{id}/sprotein/reference_based_mutations_{id}_conservation.r.ipynb",
    threads: 16
    conda:
        "../envs/R_env.yaml"
    shell:
        """
        papermill notebooks/get_reference_based_mutationsdev.ipynb \
        {log} \
        -p nextclade {input.data} \
        -p out {output.csv} \
        -p out_details {output.data} \
        """

rule analyse_contact_enrichment:
    input: 
        coservation = rez_dir + "/lineages/{id}/sprotein/{scope}_{id}_sequences_conservation_all.csv",
        contacts = config["antibody_contacts_data"] 
    output:
        data = rez_dir + "/lineages/{id}/sprotein/{scope}_{id}_contacts_vs_conservation.csv",
        image = rez_dir + "/lineages/{id}/sprotein/{scope}_{id}_contacts_review.svg" # general check view of contacts
    conda:
        "../envs/R_env.yaml"
    notebook:
        "notebooks/analyse_contacts_vs_conservation.r.ipynb"



rule get_entropies:
    input:
        #expand(rez_dir + "/lineages/{id}/sprotein/all_sequences.fasta.gz",id = ["BA.2"])
        #expand(rez_dir + "/lineages/{id}/sprotein/all_sequences_conservation_all.csv",id=['BA.2']),
        #expand(rez_dir + "/lineages/{id}/sprotein/{scope}_{id}_sequences_conservation_all.csv",id=['AY.4.5', 'Q.1', 'B.1.1.7', 'B.1.177.60', 'BA.2', 'BA.2.9'],scope=["lt","all"])
        #expand(rez_dir + "/lineages/{id}/sprotein/reference_based_mutations_{id}_data.csv",id=['AY.4.5', 'Q.1', 'B.1.1.7', 'B.1.177.60', 'BA.2', 'BA.2.9'])
        #expand(rez_dir + "/lineages/{id}/ids_lt.txt",id=['BA.2'])
        #expand(rez_dir + "/lineages/{id}/sprotein/{scope}_{id}_sequences_conservation_all.csv",id=lineages,scope=["lt","all"])
        #expand(rez_dir + "/lineages/{id}/sprotein/{scope}_{id}_contacts_vs_conservation.csv",id=lineages,scope=["lt","all"])
        data = rez_dir + "/contact_enrichment_per_position.csv"



        

