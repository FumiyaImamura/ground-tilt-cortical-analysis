#!/bin/bash
#$ -cwd
#$ -l cpu_160=1
#$ -l h_rt=1:00:00
#$ -t 1-20   # Array job with task IDs 1 to 20
#$ -o /path/to/jobs/Ridge_reductionTail_20250708/myjob.$TASK_ID.o
#$ -e /path/to/jobs/Ridge_reductionTail_20250708/myjob.$TASK_ID.e

module load matlab

# Run MATLAB without display, using a separate preferences directory
matlab -nodisplay -r "RidgeRegression_reductionTail_20250708; exit"