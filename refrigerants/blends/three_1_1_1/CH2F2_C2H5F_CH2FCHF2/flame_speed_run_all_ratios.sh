#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=2-00:00:00
#SBATCH --job-name=fcall_CH2F2_C2H5F_C2H3F3
#SBATCH --error=fc.slurm_all.log
#SBATCH --output=fc_output.slurmall.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=west


source activate cantera_env
python flame_speed_calc_all_ratios.py
