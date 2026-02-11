#!/bin/bash

#SBATCH --job-name=EXP
#SBATCH --output=EXP_%j.out
#SBATCH --error=EXP_%j.err
#SBATCH --nodelist=bigMem0
#SBATCH --time=03:00:00


module unload MATLAB
module load MATLAB/R2023b
matlab -r 'cd /home/jyliu/funm_quad; demo_exp_all_methods;'
