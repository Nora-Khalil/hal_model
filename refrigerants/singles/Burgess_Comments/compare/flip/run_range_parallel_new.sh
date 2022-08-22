#!/bin/bash
#SBATCH --job-name=flame_speeds_120-140
#SBATCH --partition=west,short
#SBATCH --exclude=c5003
#SBATCH --mem=20Gb
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=1-7%40


selected_files=(chem_1000.cti chem_0658.cti chem_0659.cti chem_0660.cti chem_0661.cti chem_0662.cti)


index=$SLURM_ARRAY_TASK_ID-1

folder_name="${selected_files[$index]}"

my_path="/work/westgroup/nora/Code/projects/halogens/refrigerants/singles/Burgess_Comments/compare/flip/flip_models/${folder_name}"


python flame_speed_calc.py $my_path
