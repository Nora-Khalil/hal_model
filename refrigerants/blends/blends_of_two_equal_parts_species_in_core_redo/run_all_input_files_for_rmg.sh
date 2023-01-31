#!/bin/sh 
#SBATCH --job-name=rmg_runs_1_21
#SBATCH --output=rmg_all.slurm.%A_%a.log
#SBATCH --error=error_rmg_all.slurm.%A_%a.log
#SBATCH --nodes=1
#SBATCH --partition=west
#SBATCH --mem=20Gb
#SBATCH --time=12-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --array=1-21

list_of_species=(CH2FCH2F_CH2F2 CH2FCH2F_CH2FCHF2 CH2FCH2F_CH3CHF2 CH2FCH2F_CH3F CH2FCHF2_CH2F2 CH2FCHF2_CH3F)

index=$SLURM_ARRAY_TASK_ID-1

folder_name="${list_of_species[$index]}"  

cd "${folder_name}"
source activate rmg_env
python-jl /home/khalil.nor/Code/RMG-Py/rmg.py input.py




