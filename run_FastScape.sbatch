#!/bin/bash -l

## General configuration options
#SBATCH -J F_TransAlp_F7
#SBATCH -e Run_4/F_Error_%j.log
#SBATCH -o Run_4/F_Screen_%j.log
#SBATCH --mail-user=paul-reinhold.eizenhoefer@uni-tuebingen.de
#SBATCH --mail-type=ALL


## Machine and CPU configuration
## Number of tasks per job:
#SBATCH -n 1
## Number of nodes:
#SBATCH -N 1
## Number of cpu's
#SBATCH -c 4

./FastScape.sh -e Run_4