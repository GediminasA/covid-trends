snakemake --use-conda --use-singularity -c 48  --configfile configs/example.yaml -k summarise_analyse_contact_enrichment 
snakemake -c 32 --use-conda --configfile configs/config.yaml --edit-notebook audines11/rez/lineages/common/common_ref_derep_data.csv
snakemake -c 32 --use-conda --configfile configs/config.yaml -f audines11/rez/pangolin_lineage_report.csv
