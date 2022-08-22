#!/bin/bash
#SBATCH --job-name=flame_speeds_120-140
#SBATCH --partition=west,short
#SBATCH --exclude=c5003
#SBATCH --mem=20Gb
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=1-21%40


# Define useful bash variables
SLURM_TASK_ID_OFFSET=0

selected_files=(copy_chem0120.cti copy_chem0121.cti copy_chem0122.cti copy_chem0123.cti copy_chem0124.cti copy_chem0125.cti copy_chem0126.cti copy_chem0127.cti copy_chem0128.cti copy_chem0129.cti copy_chem0130.cti copy_chem0131.cti copy_chem0132.cti copy_chem0133.cti copy_chem0134.cti copy_chem0135.cti copy_chem0136.cti copy_chem0137.cti copy_chem0138.cti copy_chem0139.cti copy_chem0140.cti)

index=$SLURM_ARRAY_TASK_ID-1

folder_name="${selected_files[$index]}"

my_path="/work/westgroup/nora/Code/projects/halogens/refrigerants/singles/Burgess_Comments/cantera/Nora/2_BTP/FFCM_seed/2_BTP_seed/chemkin/copies/${folder_name}"


python flame_speed_calc.py $my_path
