#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=12-00:00:00
#SBATCH --job-name={C2H5F}_mix
#SBATCH --error=tp_mix.slurm.log
#SBATCH --output=tp_mix.slurm.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=west


source activate proj_transport
PYTHONPATH=/home/khalil.nor/cantera/build/python python tp_mix_flame_speed_calc.py
##python tp_mix_flame_speed_calc.py
