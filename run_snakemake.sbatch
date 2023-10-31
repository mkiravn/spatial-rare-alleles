#!/usr/bin/env bash
#SBATCH -J simu_sfs
#SBATCH --account=pi-jnovembre
#SBATCH --mail-type ALL
#SBATCH --mail-user=mianniravn@uchicago.edu
#SBATCH --partition=jnovembre
#SBATCH --time=12:00:00

module load python
module load Anaconda3/2022.10
module load R
module load SLiM/3.0
source /software/python-anaconda-2020.02-el7-x86_64/etc/profile.d/conda.sh
source activate ./spatialrarealleles

snakemake --unlock
snakemake --profile slurm -s new_snakefile --keep-going --rerun-incomplete --printshellcmds
