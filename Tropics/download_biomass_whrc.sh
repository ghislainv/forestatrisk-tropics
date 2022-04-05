#!/bin/bash
#SBATCH --job-name=getBiomass
#SBATCH --account=agap
#SBATCH --partition=agap_normal
#SBATCH --output=/lustre/vieilledentg/getBiomass_log.txt
#SBATCH --time=24:00:00
#SBATCH --ntasks=120
#SBATCH --cpus-per-task=3

# Loading singularity module
module load singularity

# run python script inside container
singularity exec --bind /storage/replicated/cirad/projects/AMAP/vieilledentg ~/replicated/singularity_images/forestatrisk-tropics.simg python3 -u ~/Code/forestatrisk-tropics/Tropics/download_biomass_whrc.py $SLURM_JOB_ID

# End
