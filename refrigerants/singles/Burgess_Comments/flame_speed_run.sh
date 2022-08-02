#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --job-name=fc_CH4_David
#SBATCH --error=slurm/error/fc_CH4_no_halogens_slack.slurm.log
#SBATCH --output=slurm/output/fc_output_CH4_no_halogens_slack.slurm.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=short


source activate cantera_env
python flame_speed_calc.py
