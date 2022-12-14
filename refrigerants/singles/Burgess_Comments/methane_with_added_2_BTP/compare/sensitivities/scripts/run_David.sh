#!/bin/bash
#SBATCH --job-name=sens_David
#SBATCH --partition=west
#SBATCH --exclude=c5003
#SBATCH --mem=20Gb
#SBATCH --time=6-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=1



python -u sensitivity_David.py 
