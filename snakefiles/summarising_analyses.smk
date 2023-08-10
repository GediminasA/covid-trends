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
        estimates_of_stability_image = rez_dir + "/estimates_of_stability.svg"
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