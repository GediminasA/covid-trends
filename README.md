### This workflow analyses given SARS-CoV-2 data (genomic sequences, metadata of the sequences) and calculates various trends.

## Funding
This project has received funding from European Regional Development Fund (project No 13.1.1-LMT-K-718-05-0023) under grant agreement with the Research Council of Lithuania (LMTLT). Funded as European Union's measure in response to Cov-19 pandemic."

## Workflow setup

1. Install `conda`:
```bash
   wget -P miniconda https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh &&
   chmod 755 ./miniconda/Miniconda3-latest-Linux-x86_64.sh &&
   ./miniconda/Miniconda3-latest-Linux-x86_64.sh
```

2. Add path of `miniconda` to `.bashrc` if not selected option to add automatically during installation:
```bash
   cd ~/ && pth=$(pwd) &&
   echo "PATH=$PATH:$pth/miniconda3/bin" >> ~/.bashrc
```

3. Install mamba frontend for conda. This helps to handle dependancies installation.
```bash
    conda install -c conda-forge mamba
```

4. Clone the repository. nd enter it. Note - it contains submodules.
```bash 
    git clone --recurse-submodules  git@github.com:GediminasA/covid-trends.git
    cd covid-trends
```
5. Create your `conda` environment:
 ```bash
    mamba env create -f envs/covid-trends.yaml 
 ```

6. Activate created environment:
```bash
    conda activate covid-trends
```

7. Checkout exemplary data git lfs:
```bash
    git lfs fetch
    git lfs install
    git lfs checkout
```

Make sure that uidmap is install ed in the system, this could ne done like:
 ```bash
 sudo apt-get install uidmap
 ```

## Analysis run

It is best to run analysis on a computing machine which has at least:

-  250 GB RAM
-  16 cores

If small data sets are analyzed - memory requirements are lower; however, analysis of whole GISAID dataset or equivalent consumes more than 125 GB 

Analysis can be run like this:

```bash
snakemake --use-conda --use-singularity -c 48  --configfile configs/example.yaml -k assemble_raport
```
`48` - indicates number of cores to use

`configs/example.yaml` - file with input parameters

The workflow analysis with the small exemplary datasets runs ~0.5 h. The first 
run takes longer as R/julia libraries are downloaded and set up.  

## Parameter setup

Setting analysis run please use the example configuration file `configs/example.yaml` as a template and 
change as appropriate. These parameter could be changed for custom analysis:

- `work_dir` - run directory which will contains all temporary files and results 
- `sequences` - sequences in fasta format. It is recommended to analyze as much as possible, ideally, all GISAID or equivalent dataset.
- `meta` - meta data for the sequences. For acceptable format - check the example file  `test_data/four_lt_lineages_genomic_renamed_meta.csv`
- `lineages` - indicate a list for lineage-wise analysis parts, e.g "B.1.177.60	BA.2 BA.2.9 B.1.1.7"

Please leave other parameters as in the exemplary file or change them on your own risk.

