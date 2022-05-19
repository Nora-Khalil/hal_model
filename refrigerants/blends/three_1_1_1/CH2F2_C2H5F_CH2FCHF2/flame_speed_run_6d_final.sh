#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=4-00:00:00
#SBATCH --job-name=fc6_CH2F2_C2H5F_C2H3F3
#SBATCH --error=fc.slurm_6d.log
#SBATCH --output=fc_output.slurm6d.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=west


source activate cantera_env
python flame_speed_calc6d_final.py
