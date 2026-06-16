#!/bin/bash

PARTITION=$1
SETTING=$2
PARAMETERS_FILE=$3
NSEEDS=$4

module purge
module load R/4.4.0

# Make sure R can see system libraries and your personal library
export R_LIBS_USER=/apps/R/4.4.0/lib64/R/site/library:/home/acodi/Rlibs_ve_trial
export PATH=/apps/R/4.4.0/bin:/home/acodi/.local/bin:/home/acodi/bin:/home/acodi/miniconda3/bin:/home/acodi/miniconda3/condabin:/apps/bin:/usr/share/Modules/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin

#export PATH=/apps/R/4.4.0/bin:/home/acodi/.local/bin:/home/acodi/bin:/home/acodi/miniconda3/bin:/home/acodi/miniconda3/condabin:/apps/bin:/usr/share/Modules/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin

sbatch --array=1-$NSEEDS \
	--partition=$PARTITION \
	-n 1 \
	--mem-per-cpu=6G \
	--output=/projects/dbenkes/allison/shigella_vaccine_trial/scratch/${SETTING}_%A_%a.out \
	--job-name=shigella_%a \
	--export=SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID,SETTING=$SETTING,PARAMETERS_FILE=$PARAMETERS_FILE \
	--wrap "/apps/R/4.4.0/bin/Rscript R/run_sens_analysis.R"
