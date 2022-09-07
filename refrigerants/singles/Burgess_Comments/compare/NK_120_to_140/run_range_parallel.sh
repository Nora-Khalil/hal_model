#!/bin/bash
#SBATCH --job-name=one_flip_130-129
#SBATCH --partition=short
#SBATCH --exclude=c5003
#SBATCH --mem=20Gb
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=1-8%40


# Define useful bash variables
SLURM_TASK_ID_OFFSET=0

selected_files=(copy_chem_.cti copy_chem_1285.cti copy_chem_1321.cti copy_chem_1321_1285.cti copy_chem_1321_1322.cti copy_chem_1321_1322_1285.cti copy_chem_1322.cti copy_chem_1322_1285.cti)

index=$SLURM_ARRAY_TASK_ID-1

folder_name="${selected_files[$index]}"

my_path="/work/westgroup/nora/Code/projects/halogens/refrigerants/singles/Burgess_Comments/compare/handpicked/handpicked_models/copies/${folder_name}"


python flame_speed_calc.py $my_path
