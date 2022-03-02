#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ecology.ghislainv.fr
# python          :>=3
# license         :GPLv3
# ==============================================================================

# google-cloud-sdk and earthengine tool must be installed
# sudo apt-get update && sudo apt-get install google-cloud-sdk
# https://developers.google.com/earth-engine/guides/command_line

import subprocess

# Variables
# name = ["fcc_123", "prob_2020", "fcc_2050", "fcc_2100"]
name = ["fcc_2050", "fcc_2100"]
cont = ["AFR", "AME", "ASI"]
proj = "m"
proj = "aea" if proj == "a" else "merc"

# Upload from fdb server to Google Cloud Storage (GCS)
for i in name:
    for j in cont:
        filepath = ("/home/www/forestatrisk/tropics/tif/"
                    f"{i}_{j}_{proj}.tif")
        bucket = "gs://forestatrisk/tropics/"
        cmd = f"gsutil cp {filepath} {bucket}"
        subprocess.call(cmd, shell=True)

# Create image collections
for i in name:
    gee_coll = f"projects/forestatrisk/assets/v1_2020/{i}"
    cmd = f"earthengine create collection {gee_coll}"
    subprocess.call(cmd, shell=True)

# Upload from GCS to Google Earth Engine (GEE)
for i in name:
    for j in cont:
        ndval = 0 if i in ["fcc_123", "prob_2020"] else 255
        filepath = f"gs://forestatrisk/tropics/{i}_{j}_{proj}.tif"
        asset_id = f"projects/forestatrisk/assets/{i}"
        cmd = (f"earthengine upload image --asset_id={asset_id} "
               f"--pyramiding_policy=sample --nodata_value={ndval} "
               f"{filepath}")
        subprocess.call(cmd, shell=True)

# # Cleaning GEE
# for i in name:
#     cmd = ("earthengine rm -r "
#            f"projects/forestatrisk/assets/v1_2020/{i}")
#     subprocess.call(cmd, shell=True)

# EOF
