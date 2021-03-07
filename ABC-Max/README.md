# ABC-Max code
Example bash command list for setting up a project and running the pipeline to reproduce ABC-Max calculations for IBD  
Questions: Contact Jesse Engreitz at engreitz@stanford.edu

## Check out the project
    git clone git@github.com:EngreitzLab/ABC-GWAS-Paper.git
    cd ABC-GWAS-Paper/ABC-Max/

## Install and activate conda environment
    conda env create --file abc-max.yml
    conda activate abc-max

## Fetch needed data files and if needed, rename the directory
    wget --cut-dirs 4 -r ftp://ftp.broadinstitute.org/outgoing/lincRNA/Nasser2020/data/ \
    mv ftp.broadinstitute.org data

## Run snakemake pipeline
    snakemake --snakefile snakemake/ABC-Max.snakefile \
      --configfile config/ABC-Max.IBD.json \
      --config logDir=log/ outDir=out/ \
      --cores 1 \
      --keep-target-files \
      --rerun-incomplete


