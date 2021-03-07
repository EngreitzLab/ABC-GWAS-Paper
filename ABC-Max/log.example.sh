## Example bash command list for setting up a project and running the pipeline

## Check out the project
git clone git@github.com:EngreitzLab/ABC-GWAS-Paper.git
cd ABC-GWAS-Paper/ABC-Max/

## Install and activate conda environment
conda env create --file abc-max.yml
conda activate abc-max

## Fetch needed data files
wget --no-parent  -r http://mitra.stanford.edu/kundaje/projects/ABC_links/GWAS_test/Test_data/

## Run snakemake pipeline
snakemake --snakefile snakemake/ABC-Max.snakefile \
  --configfile config/ABC-Max.IBD.json \
  --config logDir=log/ outDir=out/ \
  --cores 1 \
  --keep-target-files \
  --rerun-incomplete


