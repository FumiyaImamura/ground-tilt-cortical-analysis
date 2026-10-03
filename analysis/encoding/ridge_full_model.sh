#!/bin/bash
#$ -cwd
#$ -l cpu_80=1
#$ -l h_rt=00:30:00
#$ -t 1-20   # Array job with task IDs 1 to 20
#$ -o /path/to/logs/ridge_full_model.$TASK_ID.o
#$ -e /path/to/logs/ridge_full_model.$TASK_ID.e

module load matlab

# Run MATLAB without display, using a separate preferences directory
matlab -nodisplay -r "ridge_full_model; exit"