#!/bin/bash
#SBATCH --job-name=forecast
#SBATCH --account=agap
#SBATCH --partition=agap_normal
#SBATCH --output=/lustre/vieilledentg/forecast_log.%A_%a.txt
#SBATCH --error=/lustre/vieilledentg/forecast_err.%A_%a.txt
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=3
#SBATCH --array=1-120

# Environmental variable
export MPLCONFIGDIR="/storage/replicated/cirad_users/vieilledentg/config/matplotlib"

# Loading singularity module
module load singularity

# Paths
storage="/storage/replicated/cirad/projects/AMAP/vieilledentg"
image="/storage/replicated/cirad_users/vieilledentg/singularity_images/forestatrisk-tropics.simg"
script="/home/vieilledentg/Code/forestatrisk-tropics/Tropics/forecast.py"

# run python script inside container
singularity exec --bind $storage $image python3 -u $script $SLURM_ARRAY_TASK_ID

# End
