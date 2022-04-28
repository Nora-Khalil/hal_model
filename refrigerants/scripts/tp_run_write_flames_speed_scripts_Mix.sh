#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --job-name=w_,mix
#SBATCH --error=wr_mix.slurm.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=short


PYTHONPATH=/home/khalil.nor/cantera/build/python python tp_write_flame_speed_scripts_Mix.py
