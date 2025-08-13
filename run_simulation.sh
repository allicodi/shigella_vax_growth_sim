#!/bin/bash

PARTITION=$1
SETTING=$2
PARAMETERS_FILE=$3
NSEEDS=$4

sbatch --array=1-$NSEEDS \
	--partition=$PARTITION \
	-n 1 \
	--output=/projects/dbenkes/allison/shigella_vaccine_trial/scratch/%a_%J.out \
	--job-name=shigella_%a \
	--export=SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID,SETTING=$SETTING,PARAMATERS_FALSE=$PARAMATERS_FILE \
	--wrap "/apps/R/4.4.0/bin/Rscript run_analysis.R"
