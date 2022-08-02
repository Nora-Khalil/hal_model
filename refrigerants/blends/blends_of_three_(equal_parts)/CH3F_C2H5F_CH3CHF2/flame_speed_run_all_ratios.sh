#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --job-name=Fix3fcall_CH3F(1)_C2H5F(2)_CH3CHF2(3)
#SBATCH --error=Fix3fc.slurmall.log
#SBATCH --output=Fix3fc_output.slurmall.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=short

source activate cantera_env
python flame_speed_calc_all_ratios.py

