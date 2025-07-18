#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=victor.vandermeersch@cefe.cnrs.fr
#SBATCH --job-name=gplg
#SBATCH --output=/scratch/vvandermeersch/climategrowthshifts/logleg.txt
#SBATCH --cpus-per-task=10
#SBATCH --mem=80G
#SBATCH --time 5-00:00:00

. /local/env/envconda.sh

conda activate /home/genouest/mnhn_cesco/vvandermeersch/env/rstan

Rscript run_model_legacies.R
