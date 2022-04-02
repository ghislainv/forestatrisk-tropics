#!/bin/bash
#SBATCH --job-name=rcloneSync
#SBATCH --account=agap
#SBATCH --partition=agap_normal
#SBATCH --output=/lustre/vieilledentg/rcloneSync_log.txt
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3

# Loading modules
module load singularity

# Run python script inside container
singularity exec --bind /storage/replicated/cirad/projects/AMAP/vieilledentg ~/replicated/singularity_images/forestatrisk-tropics.simg python3 -u ~/Code/forestatrisk-tropics/Backups/rclone_gdrive_meso.py

# EOF
