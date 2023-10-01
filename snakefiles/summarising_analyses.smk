rule time_trend_of_stability:
    input:
        stabilities_common = rez_dir + "/lineages/common/common_ref_derep_data.csv",
        stabilities_pairwise = rez_dir + "/pairwise_derep_results.csv",
        pangolin = rez_dir+"/pangolin_lineage_report.csv",
        meta = config["meta"],
        most_abundad_lt = rez_dir + "/most_abundand_lineages_lt.csv",
        abundance_dates_per_lineage = rez_dir + "/abundance_dates_per_lineage.csv", 
    log:
        notebook = rez_dir + "/time_trend_of_stability.ipynb"
    output:
        time_trend_of_stability_data = rez_dir + "/time_trend_of_stability.csv",
        estimates_of_stability_data = rez_dir + "/estimates_of_stability.csv",
        time_trend_of_stability_image = rez_dir + "/time_trend_of_stability.svg",
        estimates_of_stability_image = rez_dir + "/estimates_of_stability.svg",
        estimates_of_statistical_paurwise_comp_image = rez_dir + "/estimates_of_pairwise_statist.svg"
    conda:
        "../envs/R_env.yaml"
    notebook:
        "notebooks/time_trend_of_stability.r.ipynb"

rule most_abundand_lineages:
    input:
        pangolin = rez_dir+"/pangolin_lineage_report.csv",
        meta = config["meta"]
    log:
        notebook = rez_dir + "/most_abundand_lineages.ipynb"
    output:
        quarterly_most_abundand = rez_dir + "/most_abundand_lineages_quartlely_worldwide.csv",
        most_abundad_lt = rez_dir + "/most_abundand_lineages_lt.csv",
        abundance_dates_per_lineage = rez_dir + "/abundance_dates_per_lineage.csv", 
    conda:
        "../envs/R_env.yaml"
    notebook:
        "notebooks/most_abundand_lineages.r.ipynb"


rule summarise_analyse_contact_enrichment:
    input:
        lt = expand(rez_dir + "/lineages/{id}/sprotein/{scope}_{id}_contacts_vs_conservation.csv",id=lineages,scope=["lt",]),
        all = expand(rez_dir + "/lineages/{id}/sprotein/{scope}_{id}_contacts_vs_conservation.csv",id=lineages,scope=["all",]),
        abundance_dates_per_lineage = rez_dir + "/abundance_dates_per_lineage.csv", 
    log:
        notebook = rez_dir + "/contact_enrichment.ipynb"
    output:
        data = rez_dir + "/contact_enrichment.csv",
        image = rez_dir + "/contact_enrichment.svg",
    conda:
        "../envs/R_env.yaml"
    notebook:
        "../notebooks/summarise_analyse_contact_enrichment.r.ipynb"

rule summarise_analyse_contact_enrichment_per_position:
    input:
        lt = expand(rez_dir + "/lineages/{id}/sprotein/{scope}_{id}_sequences_conservation_all.csv",id=lineages,scope=["lt",]),
        all = expand(rez_dir + "/lineages/{id}/sprotein/{scope}_{id}_sequences_conservation_all.csv",id=lineages,scope=["all",]),
        abundance_dates_per_lineage = rez_dir + "/abundance_dates_per_lineage.csv", 
        contacts = config["antibody_contacts_data"] 
    output:
        data1 = rez_dir + "/contact_enrichment_per_position_1.csv",
        data2 = rez_dir + "/contact_enrichment_per_position_2.csv",
        image1 = rez_dir + "/contact_enrichment_per_position_1.svg",
        image2 = rez_dir + "/contact_enrichment_per_position_2.svg",
    params:
        lineages = lineages
    conda:
        "../envs/R_env.yaml"
    log:
        notebook = rez_dir + "/contact_enrichment_per_position.log.ipynb"
    notebook:
        "../notebooks/summarise_analyse_contact_enrichment_per_position.r.ipynb"

#fils in a template with preambe
rule assemble_raport_template:
    input:
        quarterly_most_abundand = rez_dir + "/most_abundand_lineages_quartlely_worldwide.csv",
        most_abundad_lt = rez_dir + "/most_abundand_lineages_lt.csv",
        time_trend_of_stability_data = rez_dir + "/time_trend_of_stability.csv",
        estimates_of_stability_data = rez_dir + "/estimates_of_stability.csv",
        time_trend_of_stability_image = rez_dir + "/time_trend_of_stability.svg",
        estimates_of_stability_image = rez_dir + "/estimates_of_stability.svg",
        contacts_enrichment_svg = rez_dir + "/contact_enrichment.svg",
        contacts_enrichment_csv = rez_dir + "/contact_enrichment.csv",
        estimates_of_statistical_paurwise_comp_image = rez_dir + "/estimates_of_pairwise_statist.svg",
        contact_per_pos_data1 = rez_dir + "/contact_enrichment_per_position_1.csv",
        contact_per_pos_data2 = rez_dir + "/contact_enrichment_per_position_2.csv",
        contact_per_pos_image1 = rez_dir + "/contact_enrichment_per_position_1.svg",
        contact_per_pos_image2 = rez_dir + "/contact_enrichment_per_position_2.svg",
        hap_review = rez_dir + "/info_on_cluster_and_contacts_review.png",
        hap_all_trend = rez_dir + "/info_on_cluster_and_contacts_over_time_all.svg",
        hap_lt_trend = rez_dir + "/info_on_cluster_and_contacts_over_time_lt.svg",
        hap_taus = rez_dir + "/info_on_cluster_and_contacts_over_time.csv",
        hap_haplotypes = config["work_dir"]+"/tmp/haplotypes" + "/haplotype_info_lt_per_month.csv",
        binding_antib_csv4rep = rez_dir + "/mutants_against_antibodies_rezdata.csv",
        binding_antib_png4rep = rez_dir + "/mutants_against_antibodies_rezdata.png",
        binding_ace2_csv4rep = rez_dir + "/mutants_against_ace2_rezdata.csv",
        binding_ace2_png4rep = rez_dir + "/mutants_against_ace2_rezdata.png",
    output:
        raport = tmp_dir + "/RAPORT4template.html",
    script:
        "../notebooks/raport_templategen.Rmd"


# reads preamble and creates the report
rule assemble_raport:
    input:
        raport = tmp_dir + "/RAPORT4template.html",
    output:
        raport = config["work_dir"] + "/RAPORT.html",
    script:
        "../notebooks/raport.Rmd"
