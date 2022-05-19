#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=4-00:00:00
#SBATCH --job-name=fc_CH2F2(1)_C2H5F(2)_C2H3F3(3)
#SBATCH --error=fc.slurm4_1_52.log
#SBATCH --output=fc_output.slurm4_1_52.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=west

source activate cantera_env
python flame_speed_calc4_1_52.py

