#!/bin/sh
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --job-name=flames
#SBATCH --array=0-90
#SBATCH -p express
souce activate cantera_env
python one-flame-speed.py
