#!/bin/bash

#SBATCH --partition=short
#SBATCH --mem=10G
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --output=%j_%x.log.out
#SBATCH --error=%j_%x.log.err

module load "R/4.3.2-gfbf-2023a"

Rscript --vanilla install_packages.R