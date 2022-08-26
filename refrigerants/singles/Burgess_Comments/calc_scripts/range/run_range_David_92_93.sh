#!/bin/sh

##SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --job-name=dv_92_93
##SBATCH --error=D_2BTP_130sp.error.slurm.log
##SBATCH --output=D_2BTP_130sp.output.slurm.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
#SBATCH --array=1-3
#SBATCH --partition=short


selected_files=(copy_chem0092.cti copy_chem0093.cti)

index=$SLURM_ARRAY_TASK_ID-1

file_name="${selected_files[$index]}"

my_path="/scratch/westgroup/David/halogen_models/rmg_combustion_paper_FFCM/suppressants/2-BTP/chemkin/copies/${file_name}" 
source activate cantera_env
python flame_speed_calc_range_10pts.py $my_path
