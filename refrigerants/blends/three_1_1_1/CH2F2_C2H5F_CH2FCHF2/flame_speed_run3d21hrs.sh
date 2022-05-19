#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=4-00:00:00
#SBATCH --job-name=fc_CH2F2_C2H5F_C2H3F3
#SBATCH --error=fc.slurm_3d21hrs.log
#SBATCH --output=fc_output.slurm3d21hrs.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=west


source activate cantera_env
python flame_speed_calc3d21hrs.py
