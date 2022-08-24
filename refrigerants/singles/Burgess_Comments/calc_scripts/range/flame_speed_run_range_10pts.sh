#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=4-00:00:00
#SBATCH --job-name=2_BTP_redo
#SBATCH --error=NK_2BTP_145sp.error.slurm.log
#SBATCH --output=NK_2BTP_145sp.output.slurm.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=west


source activate cantera_env
python flame_speed_calc_range_10pts_Nora.py
