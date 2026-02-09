#!/bin/bash

#SBATCH --job-name=LOG
#SBATCH --output=LOG_%j.out
#SBATCH --error=LOG_%j.err
#SBATCH --nodelist=bigMem4
#SBATCH --time=18:00:00
#SBATCH --cpus-per-task=16


module unload MATLAB
module load MATLAB/R2023b
matlab -r 'cd /home/jyliu/funm_quad; demo_log_all_methods;'
