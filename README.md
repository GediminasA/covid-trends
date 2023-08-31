mamba env create -f envs/covid-phylogeny.yam
conda activate  covid-phylogeny
snakemake -c 32 --use-conda --configfile configs/config.yaml -k  setup_R
sudo snap install julia --classic

