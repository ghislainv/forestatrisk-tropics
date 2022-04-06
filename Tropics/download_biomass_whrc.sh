#!/bin/bash
#SBATCH --job-name=getBiomass
#SBATCH --account=agap
#SBATCH --partition=agap_normal
#SBATCH --output=/lustre/vieilledentg/getBiomass_log.%A_%a.txt
#SBATCH --error=/lustre/vieilledentg/getBiomass_err.%A_%a.txt
#SBATCH --time=47:50:00
#SBATCH --cpus-per-task=3
#SBATCH --array=1-120

# Loading singularity module
module load singularity

# run python script inside container
singularity exec --bind /storage/replicated/cirad/projects/AMAP/vieilledentg ~/replicated/singularity_images/forestatrisk-tropics.simg python3 -u ~/Code/forestatrisk-tropics/Tropics/download_biomass_whrc.py $SLURM_ARRAY_TASK_ID

# End
