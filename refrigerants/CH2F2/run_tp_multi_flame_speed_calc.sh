#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=12-00:00:00
#SBATCH --job-name=mu{CH2F2}
#SBATCH --error=tp_multi.slurm.log
#SBATCH --output=tp_multi.slurm.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=west


source activate proj_transport
PYTHONPATH=/home/khalil.nor/cantera/build/python python tp_multi_flame_speed_calc.py
##python tp_multi_flame_speed_calc.py
