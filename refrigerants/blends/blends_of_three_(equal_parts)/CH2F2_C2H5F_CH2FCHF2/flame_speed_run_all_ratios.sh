#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --job-name=FIX3fcall_CH2F2_C2H5F_C2H3F3
#SBATCH --error=FIX3fc.slurm_all.log
#SBATCH --output=FIX3fc_output.slurmall.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=short


source activate cantera_env
python flame_speed_calc_all_ratios.py
