#!/bin/bash

PARTITION=$1
SETTING=$2
NSEEDS=$3

sbatch --array=1-$NSEEDS \
	--partition=$PARTITION \
	-n 1 \
	--output=/projects/dbenkes/allison/shigella/scratch/%a_%J.out \
	--job-name=shigella_%a \
	--export=SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID,SETTING=$SETTING \
	--wrap "Rscript run_simulation_cluster.R"