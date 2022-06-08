#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --job-name=test
#SBATCH --error=test.slurm.log
#SBATCH --output=test_output.slurm.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=short



python get_velocity_data.py 
