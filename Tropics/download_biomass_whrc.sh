#!/bin/bash
#SBATCH --job-name=getBiomass
#SBATCH --account=agap
#SBATCH --partition=agap_normal
#SBATCH --output=/lustre/vieilledentg/getBiomass_log.%A_%a.txt
#SBATCH --error=/lustre/vieilledentg/getBiomass_err.%A_%a.txt
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=3
#SBATCH --array=1-120

# Environmental variable
export MPLCONFIGDIR="/lustre/vieilledentg/config/matplotlib"

# Loading singularity module
module load singularity

# Paths
storage="/storage/replicated/cirad/projects/AMAP/vieilledentg"
image="/storage/replicated/cirad_users/vieilledentg/singularity_images/forestatrisk-tropics.simg"
script="/home/vieilledentg/Code/forestatrisk-tropics/Tropics/download_biomass_whrc.py"

# run python script inside container
singularity exec --bind $storage $image python3 -u $script $SLURM_ARRAY_TASK_ID

# End
