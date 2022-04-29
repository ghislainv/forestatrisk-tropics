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

# Paths
storage="/storage/replicated/cirad/projects/AMAP/vieilledentg"
image="/lustre/vieilledentg/singularity_images/forestatrisk-tropics.simg"
script="/home/vieilledentg/Code/forestatrisk-tropics/Backups/rclone_meso_gdrive.py"

# Run python script inside container
singularity exec --bind $storage $image python3 -u $script

# EOF
