#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=4-00:00:00
#SBATCH --job-name=dv_130
#SBATCH --error=D_2BTP_130sp.error.slurm.log
#SBATCH --output=D_2BTP_130sp.output.slurm.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=west


source activate cantera_env
python flame_speed_calc_range_10pts_Dv_130.py
