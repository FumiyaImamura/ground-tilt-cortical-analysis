#!/bin/bash
#$ -cwd
#$ -l cpu_40=1
#$ -l h_rt=10:00:00
#$ -t 1-20   # Array job with task IDs 1 to 20
#$ -o /path/to/jobs/TB_classi_anchor_20250712/myjob.$TASK_ID.o
#$ -e /path/to/jobs/TB_classi_anchor_20250712/myjob.$TASK_ID.e

module load matlab

# Run MATLAB without display, using a separate preferences directory
matlab -nodisplay -r "TB_classi_anchor_20250712; exit"