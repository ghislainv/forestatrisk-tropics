#!/bin/bash

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ghislainv.github.io
# license         :GPLv3
# ==============================================================================

# # gfc2020_70
# source_path="mbb:/share/nas2-amap/gvieilledent/gfc2020_70/"
# dest_path="/home/forestatrisk-tropics/gfc2020_70/"
# data_raw="*/data_raw/*"
# rsync -e "ssh -i /home/ghislain/.ssh/id_rsa" -av --exclude $data_raw --delete-after --stats --progress $source_path $dest_path

# jrc2020
source_path="meso:/storage/replicated/cirad/projects/AMAP/vieilledentg/jrc2020/"
dest_path="/home/forestatrisk-tropics/jrc2020/"
data_raw="*/data_raw/*"
rsync -e "ssh -i /home/ghislain/.ssh/id_rsa_ghvi" -av --exclude $data_raw --delete-after --stats --progress $source_path $dest_path

# home meso
source_path="meso:/home/vieilledentg/"
dest_path="/home/forestatrisk-tropics/home_meso/"
rsync -e "ssh -i /home/ghislain/.ssh/id_rsa_ghvi" -av --exclude=".singularity" --exclude=".cache" --delete-after --stats --progress $source_path $dest_path

# EOF
