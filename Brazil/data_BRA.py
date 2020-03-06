#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ecology.ghislainv.fr
# python_version  :>=2.7
# license         :GPLv3
# ==============================================================================

# Imports
import os
from shutil import copy2  # To copy files
import forestatrisk as far
import subprocess
import pywdpa
from dotenv import load_dotenv
load_dotenv("/home/ghislain/Code/forestatrisk-tropics/.env")
import ee
ee.Initialize()

# Isocode for Brazil states
iso3 = ["BRA-AC", "BRA-AL", "BRA-AM", "BRA-AP", "BRA-BA", "BRA-CE", "BRA-ES",
        "BRA-GO", "BRA-MA", "BRA-MG", "BRA-MS", "BRA-MT", "BRA-PA", "BRA-PB",
        "BRA-PE", "BRA-PI", "BRA-PR", "BRA-RJ", "BRA-RN", "BRA-RO", "BRA-RR",
        "BRA-RS", "BRA-SC", "BRA-SE", "BRA-SP", "BRA-TO"]
nstates = len(iso3)

# Dissolve
output_f = "_gadm36_BRA_1.shp"
cmd = "ogr2ogr -f 'ESRI Shapefile' -dialect SQLITE \
      -sql 'SELECT ST_union(ST_buffer(Geometry,0.001)),GID_1 \
      FROM gadm36_BRA_2 GROUP BY GID_1' " + output_f + " gadm36_BRA_2.shp"
subprocess.call(cmd, shell=True)

# Merge distrito-federal (BRA.9_1) and goias (BRA.7_1)
output_f = "gadm36_BRA_1.shp"
cmd = "ogr2ogr -f 'ESRI Shapefile' -dialect SQLITE \
      -sql 'SELECT ST_union(ST_buffer(Geometry,0.001)),GID_1 \
      FROM _gadm36_BRA_1 GROUP BY GID_1' " + output_f + " _gadm36_BRA_1.shp"
subprocess.call(cmd, shell=True)

# Brazil PA
pywdpa.get_token()
pywdpa.get_wdpa("AUS")

# State boundaries
for i in range(nstates):
    iso = iso3[i]
    far.make_dir(iso + "/data_raw")
    output_f = iso + "/data_raw/" + "gadm36_" + iso + "_0.shp"
    cmd = "ogr2ogr -f 'ESRI Shapefile' -dialect SQLITE \
    -sql 'SELECT * FROM gadm36_BRA_1 WHERE State_ID=" + '"' + iso + '"' + "' " + output_f + " gadm36_BRA_1.shp"
    subprocess.call(cmd, shell=True)

# Copy Brasil PA to each state folder
for i in range(nstates):
    iso = iso3[i]
    copy2("pa_BRA.dbf", iso + "/data_raw" + "/pa_" + iso + ".dbf")
    copy2("pa_BRA.prj", iso + "/data_raw" + "/pa_" + iso + ".prj")
    copy2("pa_BRA.shp", iso + "/data_raw" + "/pa_" + iso + ".shp")
    copy2("pa_BRA.shx", iso + "/data_raw" + "/pa_" + iso + ".shx")

# Run forest computation
for i in range(nstates):
    iso = iso3[i]
    #shp_name = iso + "/data_raw/" + "gadm36_" + iso + "_0.shp"
    #far.data.extent_shp(shp_name)
    far.data.country_forest_run(
        iso,
        proj="EPSG:3395",
        output_dir=iso + "/data_raw",
        keep_dir=True,
        fcc_source="jrc",
        perc=50,
        gdrive_remote_rclone="gdrive_gv",
        gdrive_folder="GEE-forestatrisk-tropics")

# EOF
