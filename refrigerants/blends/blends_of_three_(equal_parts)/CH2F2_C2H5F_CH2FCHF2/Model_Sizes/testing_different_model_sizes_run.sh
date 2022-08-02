#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --job-name=diff_models
#SBATCH --error=testingdifferentmodels.slurm.log
#SBATCH --output=testingdifferentmodels.outputslurm.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=short


source activate cantera_env 
python testing_different_model_sizes.py
