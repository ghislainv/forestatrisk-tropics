#!/bin/bash
#SBATCH --job-name=FARforecast
#SBATCH --account=agap
#SBATCH --partition=agap_normal
#SBATCH --output=/lustre/vieilledentg/getBiomass_log.%A_%a.txt
#SBATCH --error=/lustre/vieilledentg/getBiomass_err.%A_%a.txt
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=3
#SBATCH --array=1-120

# Environmental variable
export MPLCONFIGDIR="/storage/replicated/cirad_users/vieilledentg/config/matplotlib"

# Loading singularity module
module load singularity

# run python script inside container
singularity exec --bind /storage/replicated/cirad/projects/AMAP/vieilledentg ~/replicated/singularity_images/forestatrisk-tropics.simg python3 -u ~/Code/forestatrisk-tropics/Tropics/forecast.py $SLURM_ARRAY_TASK_ID

# End
