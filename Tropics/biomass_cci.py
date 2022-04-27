#!/usr/bin/python
# -*- coding: utf-8 -*-

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ghislainv.github.io
# python_version  :2.7
# license         :GPLv3
# ==============================================================================

# PRIOR TO EXECUTING THE FOLLOWING SCRIPT, YOU MUST HAVE CREDENTIALS FOR
# 1. Google Earth Engine
# 2. rclone with Google Drive: https://rclone.org/drive/
# 3. WDPA: https://www.protectedplanet.net/

# Standard library imports
import os
import pkg_resources
import shutil  # for rmtree
import sys
import subprocess

# Third party imports
import forestatrisk as far
import pandas as pd
from osgeo import gdal

# Arguments
index_ctry = int(sys.argv[1]) - 1

os.environ["PROJ_LIB"] = "/home/ghislain/.pyenv/versions/miniconda3-latest/envs/conda-far/share/proj"

# ==================
# Settings
# GDAL
os.environ["GDAL_CACHEMAX"] = "1024"
# ==================

# Directories
cluster = "meso"  # (can be fdb, mbb, meso)
if cluster == "meso":
    home_dir = "/home/vieilledentg/"
    work_dir = "/storage/replicated/cirad/projects/AMAP/vieilledentg/jrc2020/"
    temp_dir = "/lustre/vieilledentg/tmp/"
if cluster == "mbb":
    home_dir = "/home/gvieilledent/"
    work_dir = "/share/nas2-amap/gvieilledent/jrc2020/"
    temp_dir = "/share/nas2-amap/gvieilledent/tmp/"
if cluster == "fdb":
    home_dir = "/home/ghislain/"
    work_dir = "/home/forestatrisk-tropics/jrc2020/"
    temp_dir = "/home/forestatrisk-tropics/tmp/"

# Country isocode
file_ctry_run = pkg_resources.resource_filename("forestatrisk",
                                                "data/ctry_run.csv")
data_ctry_run = pd.read_csv(file_ctry_run, sep=";", header=0)
iso3 = list(data_ctry_run.iso3)
nctry = len(iso3)  # 120


# progress bar
def progress_callback(complete, message, callback_data):
    '''Emit progress report in numbers for 10% intervals and dots for 3%'''
    if int(complete*100) % 10 == 0:
        print(f'{complete*100:.0f}', end='', flush=True)
    elif int(complete*100) % 3 == 0:
        print(f'{callback_data}', end='', flush=True)


# Function for multiprocessing
def run_country(iso3):

    # GDAL temp directory
    far.make_dir(temp_dir + "tmp_" + iso3)
    os.environ["CPL_TMPDIR"] = temp_dir + "tmp_" + iso3

    # Set original working directory
    cont = data_ctry_run.cont_run[data_ctry_run["iso3"] == iso3].iloc[0]
    owd = work_dir + cont
    os.chdir(owd)
    far.make_dir(iso3)
    os.chdir(os.path.join(owd, iso3))

    # Albers Equal Area projections
    if cont == "Africa":
        proj = ("+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 "
                "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
    elif cont == "Asia":
        proj = ("+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 "
                "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
    else:
        proj = ("+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 "
                "+x_0=0 +y_0=0 +ellps=aust_SA +units=m no_defs")

    # Extent
    extent_proj = far.data.extent_shp("data/ctry_PROJ.shp")
    xmin_reg = np.floor(extent_proj[0] - 5000)
    ymin_reg = np.floor(extent_proj[1] - 5000)
    xmax_reg = np.ceil(extent_proj[2] + 5000)
    ymax_reg = np.ceil(extent_proj[3] + 5000)
    extent_reg = [xmin_reg, ymin_reg, xmax_reg, ymax_reg]
    # extent_reg_str = " ".join(map(str, extent_reg))

    # Resampling CCI data without compression using .vrt file
    input_file = home_dir + ("Code/forestatrisk-tropics/AGB/"
                             "CCI_v3/biomass_cci.tif")
    output_file = "data/emissions/biomass_cci_wrap.vrt"
    # gdal_cmd = ["gdalwarp",
    #             "-overwrite -tap -multi",
    #             "-tr 100 100",
    #             "-r bilinear",
    #             "-of VRT",
    #             "-wo NUM_THREADS=ALL_CPUS",
    #             "-te", extent_reg_str,
    #             "-t_srs", "'" + proj + "'",
    #             input_file, output_file]
    # subprocess.call(" ".join(gdal_cmd), shell=True)

    param = gdal.WarpOptions(options=["overwrite", "tap"],
                             format="VRT",
                             xRes=100, yRes=100,
                             outputBounds=extent,
                             srcSRS="EPSG:4326", dstSRS=proj,
                             resampleAlg=gdal.GRA_Bilinear,
                             multithread=True,
                             warpOptions=["NUM_THREADS=ALL_CPUS"])
    gdal.Warp(output_file, input_file, options=param)

    # Compressing
    input_file = "data/emissions/biomass_cci_wrap.vrt"
    output_file = "data/emissions/biomass_cci.tif"
    param = gdal.TranslateOptions(options=["overwrite"],
                                  format="GTiff",
                                  creationOptions=["COMPRESS=LZW",
                                                   "PREDICTOR=2",
                                                   "BIGTIFF=YES"],
                                  callback=progress_callback,
                                  callback_data=".")
    gdal.Translate(output_file, input_file, options=param)

    # Remove GDAL temp directory
    shutil.rmtree(temp_dir + "tmp_" + iso3)

    # Print country iso code
    print(iso3)


# Run country
# for i in range(nctry):
#     run_country(iso3[i])

run_country(iso3[index_ctry])

# End
