#!/bin/bash
#$ -cwd
#$ -l node_q=1
#$ -l h_rt=1:30:00
#$ -t 1-20   # Array job with task IDs 1 to 20
#$ -o /path/to/jobs/Ridge_cont_eigpos_20250628/myjob.$TASK_ID.o
#$ -e /path/to/jobs/Ridge_cont_eigpos_20250628/myjob.$TASK_ID.e

module load matlab

# Run MATLAB without display, using a separate preferences directory
matlab -nodisplay -r "RidgeRegression_eigposture_20250628; exit"