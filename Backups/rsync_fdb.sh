#!/bin/bash

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ghislainv.github.io
# license         :GPLv3
# ==============================================================================

source_path="mbb:/share/nas2-amap/gvieilledent/gfc2020_70/"
dest_path="/home/forestatrisk-tropics/gfc2020_70/"
data_raw="*/data_raw/*"
rsync -e "ssh -i /home/ghislain/.ssh/id_rsa" -av --exclude $data_raw --delete-after --stats --progress $source_path $dest_path

# EOF
