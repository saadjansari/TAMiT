#!/bin/bash

#BATCH --job-name=simCell
#SBATCH --qos=condo
#SBATCH --partition=shas
#SBATCH --account=ucb-summit-smr
#SBATCH --output=fit.log
#SBATCH --error=fit.err
#SBATCH --time=44:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12

export SCRATCH=/scratch/summit/saan8193
mkdir -p $SCRATCH/$SLURM_JOB_ID

module purge
module load matlab/R2018b

matlab -nodesktop -r "clear all; setenv('TZ','America/Denver'); cd '/projects/saan8193/ImageAnalysis/SingleCell'; simAndFitBatch('/scratch/summit/saan8193/SingleCell/simCell_1')"

rm -rf $SCRATCH/$SLURM_JOB_ID
