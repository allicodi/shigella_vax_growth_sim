#!/bin/bash

PARTITION=$1
SETTING=$2
TRUTH_DIR=$3

module purge
module load R/4.4.0

# Make sure R can see system libraries and your personal library
export R_LIBS_USER=/apps/R/4.4.0/lib64/R/site/library:/home/acodi/Rlibs_ve_trial
export PATH=/apps/R/4.4.0/bin:/home/acodi/.local/bin:/home/acodi/bin:/home/acodi/miniconda3/bin:/home/acodi/miniconda3/condabin:/apps/bin:/usr/share/Modules/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin

sbatch \
    --partition=$PARTITION \
    -n 1 \
    --mem-per-cpu=6G \
    --output=/projects/dbenkes/allison/shigella_vaccine_trial/scratch/truth_${SETTING}_%J.out \
    --job-name=truth_${SETTING} \
    --wrap "/apps/R/4.4.0/bin/Rscript R/get_sens_truth.R $SETTING $TRUTH_DIR"