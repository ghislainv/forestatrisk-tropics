#!/bin/bash
#SBATCH --job-name=build_simg
#SBATCH --account=agap
#SBATCH --partition=agap_normal
#SBATCH --output=/lustre/vieilledentg/build_simg_log.txt
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=12Go

# Environment variables
export SINGULARITY_TMPDIR=/home/vieilledentg/scratch
export SINGULARITY_CACHEDIR=/home/vieilledentg/scratch
unset XDG_RUNTIME_DIR

# Loading modules
module load singularity
module load squashfs

# Build from docker image
singularity build ~/scratch/forestatrisk-tropics.simg docker://ghislainv/docker-forestatrisk-tropics

# Clean scratch
mv ~/scratch/forestatrisk-tropics.simg ~/replicated/singularity_images/

# End
