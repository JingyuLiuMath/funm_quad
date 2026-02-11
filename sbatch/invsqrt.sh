#!/bin/bash

#SBATCH --job-name=INVSQRT
#SBATCH --output=INVSQRT_%j.out
#SBATCH --error=INVSQRT_%j.err
#SBATCH --nodelist=bigMem0
#SBATCH --time=03:00:00


module unload MATLAB
module load MATLAB/R2023b
matlab -r 'cd /home/jyliu/funm_quad; demo_invsqrt_all_methods;'
