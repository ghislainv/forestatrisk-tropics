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
data="*/data/*"
rsync -e ssh -av --exclude {$data_raw, $data} --delete-after --stats --progress --dry-run $source_path $dest_path

# EOF
