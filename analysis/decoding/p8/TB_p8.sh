#!/bin/bash
#$ -cwd
#$ -l cpu_80=1
#$ -l h_rt=20:00:00
#$ -t 1-20   # Array job with task IDs 1 to 20
#$ -o /path/to/jobs/TBreg_p8_20250720/myjob.$TASK_ID.o
#$ -e /path/to/jobs/TBreg_p8_20250720/myjob.$TASK_ID.e

module load matlab

# Run MATLAB without display, using a separate preferences directory
matlab -nodisplay -r "TB_regress_p8_20250719; exit"