#!/bin/bash
#$ -cwd
#$ -l cpu_40=1
#$ -l h_rt=24:00:00
#$ -t 1-20   # Array job with task IDs 1 to 20
#$ -o /path/to/jobs/TBregress_multi_allreg_anchor_20250617_CSP2/myjob.$TASK_ID.o
#$ -e /path/to/jobs/TBregress_multi_allreg_anchor_20250617_CSP2/myjob.$TASK_ID.e

module load matlab

# Run MATLAB without display, using a separate preferences directory
matlab -nodisplay -r "TB_regress_multireg_anchor_20250617; exit"