snakemake -c 32 --use-conda --configfile configs/config.yaml -f audines11/rez/pangolin_lineage_report.csv
snakemake --use-conda --use-singularity -c 48  --configfile configs/example.yaml -k assemble_raport
