#!/bin/sh
#set a job name
#SBATCH --job-name=flame_rmg.%a

#a file for job output, you can check job progress
#SBATCH --output=%a.log

# a file for errors from the job
#SBATCH --error=%a.slurm.log

#time you think you need; default is one day
# d-hh:mm:ss
#SBATCH --time=24:00:00
#SBATCH --exclude=c5003
#number of tasks you are requesting
#SBATCH -N 1
#SBATCH -n 6
#SBATCH --mem=30GB
##SBATCH --ntasks-per-node=1
##SBATCH --exclusive

#partition to use
#SBATCH --partition=short,west

#SBATCH --array=1-200%100
##SBATCH --array=1
##SBATCH --array=24,39,27,69,42,7,94,57,33,17,97,20,63,28,99,46,13,23,35,67,6,22,10,26,80,16,1,37,8,87,79,11,88,47,75,36,77

source activate rmg_env
export AUTOTST=../AutoTST
export halogen_data=../halogen_data
python flame_array.py
