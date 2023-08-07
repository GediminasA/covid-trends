rule time_trend_of_stability:
    input:
        stabilities_common = rez_dir + "/lineages/common/common_ref_derep_data.csv",
        stabilities_pairwise = rez_dir + "/pairwise_derep_results.csv",
        pangolin = rez_dir+"/pangolin_lineage_report.csv",
        meta = config["meta"]
    log:
        notebook = rez_dir + "/time_trend_of_stability.ipynb"
    output:
        stabilities_pairwise = rez_dir + "/time_trend_of_stability.csv"
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
        stabilities_pairwise = rez_dir + "/most_abundand_lineages.csv"
    conda:
        "../envs/R_env.yaml"
    notebook:
        "notebooks/most_abundand_lineages.r.ipynb"