#!/bin/bash
#SBATCH --job-name=fc_2blends
#SBATCH --output=fc_all.slurm.%x.log
#SBATCH --error=error_fc_all.slurm.%x.log
#SBATCH --nodes=1
#SBATCH --partition=short
#SBATCH --mem=20Gb
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --array=1-21

list_of_species_folders=(CH2FCHF2_CH3F C2H5F_CH2F2 C2H5F_CH3F CH2FCH2F_CH3F CH3CF3_CH2FCH2F CH2F2_CH3F CH3CF3_CH2FCHF2 CH2F2_CH3CF3 CH3CHF2_CH2FCHF2 CH2FCH2F_CH2FCHF2 CH3CHF2_CH2F2 CH3CF3_CH3CHF2 CH3F_CH3CHF2 CH3F_CH3CF3 CH2FCH2F_CH2F2 C2H5F_CH3CHF2 CH2FCH2F_CH3CHF2 C2H5F_CH2FCHF2 C2H5F_CH2FCH2F CH2FCHF2_CH2F2 CH3CF3_C2H5F)

#list_of_species_folders=(CH3F_CH3CF3)

index=$SLURM_ARRAY_TASK_ID-1

folder_name="${list_of_species_folders[$index]}"  

my_path="/work/westgroup/nora/Code/projects/halogens/refrigerants/blends/${folder_name}/chemkin"

echo $index
echo $my_path
cd $my_path
source activate cantera_env
python converter.py