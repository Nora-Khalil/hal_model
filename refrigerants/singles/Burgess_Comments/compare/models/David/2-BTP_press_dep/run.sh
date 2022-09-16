#!/bin/sh

#number of tasks you are requesting
#SBATCH --nodes=1
##SBATCH --exclude=c3040
#SBATCH --job-name=2-BTP
##SBATCH --exclusive
#SBATCH -o rmg.out
#SBATCH --error=rmg.slurm.log
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=8Gb
#SBATCH --ntasks=1 

#partition to use.
#SBATCH --partition=west

#an array for the job.
#SBATCH --array=1

####################################################
source activate rmg_env_jl
export RMGpy=/scratch/westgroup/David/halogen_models/rmg/RMG-Py
export PYTHONPATH=$RMGpy:$PYTHONPATH
srun python-jl $RMGpy/rmg.py -n 5 input.py
