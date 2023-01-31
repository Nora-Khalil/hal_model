#!/bin/sh 
#SBATCH --job-name=rmg_runs_1_21
#SBATCH --output=rmg_all.slurm.%A_%a.log
#SBATCH --error=error_rmg_all.slurm.%A_%a.log
#SBATCH --nodes=1
#SBATCH --partition=short
#SBATCH --mem=20Gb
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --array=1-21

list_of_species=(C2H5F_CH2F2 C2H5F_CH2FCH2F C2H5F_CH2FCHF2 C2H5F_CH3CHF2 C2H5F_CH3F CH2F2_CH3CF3 CH2F2_CH3F CH2FCH2F_CH2F2 CH2FCH2F_CH2FCHF2 CH2FCH2F_CH3CHF2 CH2FCH2F_CH3F CH2FCHF2_CH2F2 CH2FCHF2_CH3F CH3CF3_C2H5F CH3CF3_CH2FCH2F CH3CF3_CH2FCHF2 CH3CF3_CH3CHF2 CH3CHF2_CH2F2 CH3CHF2_CH2FCHF2 CH3F_CH3CF3 CH3F_CH3CHF2)

index=$SLURM_ARRAY_TASK_ID-1

folder_name="${list_of_species[$index]}"  

cd "${folder_name}/chemkin"


converter="/work/westgroup/nora/Code/projects/halogens/refrigerants/blends/blends_of_two_equal_parts_species_in_core_redo/${folder_name}/chemkin/converter_2nd_try.py"
source activate cantera_env
python $converter



