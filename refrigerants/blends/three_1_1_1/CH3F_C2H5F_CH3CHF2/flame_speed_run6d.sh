#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=3-00:00:00
#SBATCH --job-name=fc_CH3F(1)_C2H5F(2)_CH3CHF2(3)
#SBATCH --error=fc.slurm6d.log
#SBATCH --output=fc_output.slurm6d.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=west

source activate cantera_env
python flame_speed_calc6d.py

