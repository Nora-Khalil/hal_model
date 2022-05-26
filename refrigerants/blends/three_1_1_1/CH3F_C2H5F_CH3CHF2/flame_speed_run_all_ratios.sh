#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=2-00:00:00
#SBATCH --job-name=fcall_CH3F(1)_C2H5F(2)_CH3CHF2(3)
#SBATCH --error=fc.slurmall.log
#SBATCH --output=fc_output.slurmall.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=west

source activate cantera_env
python flame_speed_calc_all_ratios.py

