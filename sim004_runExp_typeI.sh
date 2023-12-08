#!/bin/bash
# Job name:
#SBATCH --job-name=FEMA_TypeI
#
# Project:
#SBATCH --account=p697 --partition=bigmem
#
# Wall clock limit:
#SBATCH --time=150:00:00
#
# Max memory usage:
#SBATCH --mem-per-cpu=20G
#
# Number of cores:
#SBATCH --cpus-per-task=40

## Set up job environment:
module purge   # clear any inherited modules
module load MATLAB/2023a
set -o errexit # exit on errors

## Do some work:
cd /ess/p697/cluster/users/parekh/2023-02-02_FEMA-Experiments/2023-11-17_Redone/scripts
matlab -nodisplay -nosplash -nodesktop < sim004_doExp_typeI.m