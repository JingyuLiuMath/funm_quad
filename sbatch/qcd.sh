#!/bin/bash

#SBATCH --job-name=qcd
#SBATCH --output=qcd_%j.out
#SBATCH --error=qcd_%j.err
#SBATCH --nodelist=bigMem1
#SBATCH --cpus-per-task=64

module unload MATLAB
module load MATLAB/R2023b

echo "=========================================="
echo "        SLURM JOB INFORMATION"
echo "=========================================="
echo "Job ID:           $SLURM_JOB_ID"
echo "Job Name:         $SLURM_JOB_NAME"
echo "Submit Directory: $SLURM_SUBMIT_DIR"
echo "Submit Host:      $SLURM_SUBMIT_HOST"
echo "Partition:        $SLURM_JOB_PARTITION"
echo "Node List:        $SLURM_JOB_NODELIST"
echo "Nodes:            $SLURM_NNODES"
echo "Tasks:            $SLURM_NTASKS"
echo "CPUs per Task:    $SLURM_CPUS_PER_TASK"
echo "Total CPUs:       $SLURM_NPROCS"
echo "Working Directory: $(pwd)"
echo "=========================================="
echo ""
echo "        RESOURCE ALLOCATION"
echo "=========================================="
echo "Memory per Node:  $SLURM_MEM_PER_NODE"
echo "Time Limit:       $SLURM_TIMELIMIT"
echo "=========================================="
echo ""
echo "        ENVIRONMENT VARIABLES"
echo "=========================================="
echo "PATH:             $PATH"
echo "LD_LIBRARY_PATH:  $LD_LIBRARY_PATH"
echo "LIBRARY_PATH:     $LIBRARY_PATH"
echo "CPATH:            $CPATH"
echo "=========================================="
echo ""
echo "        MODULE INFORMATION"
echo "=========================================="
module list 2>&1
echo "=========================================="
echo ""
echo "Job started at: $(date)"
echo "=========================================="
echo ""

matlab -r 'cd /home/jyliu/funm_quad; test_qcd;'

echo ""
echo "MATLAB finished at: $(date)"
echo "=========================================="
echo "Job completed"
echo "=========================================="
