#!/bin/bash
#SBATCH --job-name=rmg_runs_1_42
#SBATCH --output=rmg_all.slurm.%A_%a.log
#SBATCH --error=error_rmg_all.slurm.%A_%a.log
#SBATCH --nodes=1
#SBATCH --partition=west
#SBATCH --mem=20Gb
#SBATCH --time=6-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --array=1-4

##list_of_species = (CH2FCHF2_CH3F C2H5F_CH2F2 C2H5F_CH3F CH2FCH2F_CH3F CH3CF3_CH2FCH2F CH2F2_CH3F CH3CF3_CH2FCHF2 CH2F2_CH3CF3 CH3CHF2_CH2FCHF2 CH2FCH2F_CH2FCHF2 CH3CHF2_CH2F2 CH3CF3_CH3CHF2 CH3F_CH3CHF2 CH3F_CH3CF3 CH2FCH2F_CH2F2 C2H5F_CH3CHF2 CH2FCH2F_CH3CHF2 C2H5F_CH2FCHF2 C2H5F_CH2FCH2F CH2FCHF2_CH2F2 CH3CF3_C2H5F)

list_of_species_folders=(CH2FCHF2_CH3F C2H5F_CH2F2 C2H5F_CH3F)

index=$SLURM_ARRAY_TASK_ID-1

folder_name="${list_of_species_folders[$index]}"  


cd "${folder_name}"
source activate rmg_env
python-jl /home/khalil.nor/Code/RMG-Py/rmg.py input.py




