#!/bin/bash
#SBATCH --job-name=one_flip_130-129
#SBATCH --partition=west
#SBATCH --exclude=c5003
#SBATCH --mem=20Gb
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=1


# Define useful bash variables
SLURM_TASK_ID_OFFSET=0

selected_files=(chem_869_840.cti)

index=$SLURM_ARRAY_TASK_ID-1

folder_name="${selected_files[$index]}"

my_path="/work/westgroup/nora/Code/projects/halogens/refrigerants/singles/Burgess_Comments/compare/handpicked/handpicked_models/${folder_name}"


python -u  flame_speed_calc-Copy1.py $my_path


