snakemake -c 32 --use-conda --configfile configs/config.yaml --edit-notebook audines11/rez/lineages/common/common_ref_derep_data.csv
snakemake -c 32 --use-conda --configfile configs/config.yaml -f audines11/rez/pangolin_lineage_report.csv
