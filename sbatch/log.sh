#!/bin/bash

#SBATCH --job-name=LOG
#SBATCH --output=LOG_%j.out
#SBATCH --error=LOG_%j.err
#SBATCH --nodelist=bigMem0
#SBATCH --time=00:30:00


module unload MATLAB
module load MATLAB/R2023b
matlab -r 'cd /home/jyliu/funm_quad; demo_log_all_methods;'
