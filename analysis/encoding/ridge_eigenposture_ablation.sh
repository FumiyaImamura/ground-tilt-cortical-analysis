#!/bin/bash
#$ -cwd
#$ -l node_q=1
#$ -l h_rt=1:30:00
#$ -t 1-20   # Array job with task IDs 1 to 20
#$ -o /path/to/logs/ridge_eigenposture_ablation.$TASK_ID.o
#$ -e /path/to/logs/ridge_eigenposture_ablation.$TASK_ID.e

module load matlab

# Run MATLAB without display, using a separate preferences directory
matlab -nodisplay -r "ridge_eigenposture_ablation; exit"