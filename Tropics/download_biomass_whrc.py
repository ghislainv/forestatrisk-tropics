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

# Standard library imports
import os
import pkg_resources
# import shutil  # for rmtree
# import sys

# Third party imports
import ee
import forestatrisk as far
import pandas as pd

# index_ctry = int(sys.argv[1])-1

# ==================
# Settings
# Earth engine
ee.Initialize()
# GDAL
os.environ["GDAL_CACHEMAX"] = "1024"
# ==================

# Country isocode
file_ctry_run = pkg_resources.resource_filename("forestatrisk",
                                                "data/ctry_run.csv")
data_ctry_run = pd.read_csv(file_ctry_run, sep=";", header=0)
iso3 = list(data_ctry_run.iso3)
nctry = len(iso3)  # 120

# Cluster (can be fdb, mbb, meso)
cluster = "meso"

# Directories
if cluster == "meso":
    work_dir = "~/projects/AMAP/vieilledentg/jrc2020/"
    temp_dir = "~/scratch/tmp/"
if cluster == "mbb":
    work_dir = "/share/nas2-amap/gvieilledent/jrc2020/"
    temp_dir = "/share/nas2-amap/gvieilledent/tmp/"
if cluster == "fdb":
    work_dir = "/home/forestatrisk-tropics/jrc2020/"

# Projections
proj_aea_Afr = ("+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 "
                "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
proj_aea_Asi = ("+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 "
                "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
proj_aea_Ame = ("+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 "
                "+x_0=0 +y_0=0 +ellps=aust_SA +units=m no_defs")
proj_aea_Afr_wkt = ('PROJCS["Africa_Albers_Equal_Area_Conic",'
                    'GEOGCS["GCS_WGS_1984",'
                    'DATUM["WGS_1984",'
                    'SPHEROID["WGS_1984",6378137,298.257223563]],'
                    'PRIMEM["Greenwich",0],'
                    'UNIT["Degree",0.017453292519943295]],'
                    'PROJECTION["Albers_Conic_Equal_Area"],'
                    'PARAMETER["False_Easting",0],'
                    'PARAMETER["False_Northing",0],'
                    'PARAMETER["longitude_of_center",25],'
                    'PARAMETER["Standard_Parallel_1",20],'
                    'PARAMETER["Standard_Parallel_2",-23],'
                    'PARAMETER["latitude_of_center",0],'
                    'UNIT["Meter",1],'
                    'AUTHORITY["EPSG","102022"]]')
proj_aea_Asi_wkt = ('PROJCS["Asia_South_Albers_Equal_Area_Conic",'
                    'GEOGCS["GCS_WGS_1984",'
                    'DATUM["WGS_1984",'
                    'SPHEROID["WGS_1984",6378137,298.257223563]],'
                    'PRIMEM["Greenwich",0],'
                    'UNIT["Degree",0.017453292519943295]],'
                    'PROJECTION["Albers_Conic_Equal_Area"],'
                    'PARAMETER["False_Easting",0],'
                    'PARAMETER["False_Northing",0],'
                    'PARAMETER["longitude_of_center",125],'
                    'PARAMETER["Standard_Parallel_1",7],'
                    'PARAMETER["Standard_Parallel_2",-32],'
                    'PARAMETER["latitude_of_center",-15],'
                    'UNIT["Meter",1],'
                    'AUTHORITY["EPSG","102028"]]')
proj_aea_Ame_wkt = ('PROJCS["South_America_Albers_Equal_Area_Conic",'
                    'GEOGCS["GCS_South_American_1969",'
                    'DATUM["South_American_Datum_1969",'
                    'SPHEROID["GRS_1967_Truncated",6378160,298.25]],'
                    'PRIMEM["Greenwich",0],'
                    'UNIT["Degree",0.017453292519943295]],'
                    'PROJECTION["Albers_Conic_Equal_Area"],'
                    'PARAMETER["False_Easting",0],'
                    'PARAMETER["False_Northing",0],'
                    'PARAMETER["longitude_of_center",-60],'
                    'PARAMETER["Standard_Parallel_1",-5],'
                    'PARAMETER["Standard_Parallel_2",-42],'
                    'PARAMETER["latitude_of_center",-32],'
                    'UNIT["Meter",1],'
                    'AUTHORITY["EPSG","102033"]]')


# Function for multiprocessing
def run_country(iso3):

    # # GDAL temp directory
    # far.make_dir(temp_dir + "tmp_" + iso3)
    # os.environ["CPL_TMPDIR"] = temp_dir + "tmp_" + iso3

    # Set original working directory
    cont = data_ctry_run.cont_run[data_ctry_run["iso3"] == iso3].iloc[0]
    owd = work_dir + cont
    os.chdir(owd)
    far.make_dir(iso3)
    os.chdir(os.path.join(owd, iso3))

    # Albers Equal Area projections
    if cont == "Africa":
        proj_aea = proj_aea_Afr
        proj_wkt = proj_aea_Afr_wkt
    elif cont == "Asia":
        proj_aea = proj_aea_Asi
        proj_wkt = proj_aea_Afr_wkt
    else:
        proj_aea = proj_aea_Ame
        proj_wkt = proj_aea_Afr_wkt

    # Extract data for country on GEE
    far.data.country_biomass_run(
        iso3=iso3,
        proj=proj_wkt,
        output_dir="data_raw",
        keep_dir=True,
        gdrive_remote_rclone="gdrive_gv",
        gdrive_folder="GEE_biomass_whrc")

    # # Download data locally from Google Drive
    # far.data.country_biomass_download(
    #     iso3=iso3,
    #     gdrive_remote_rclone="gdrive_gv",
    #     gdrive_folder="GEE_biomass_whrc",
    #     output_dir="data_raw")

    # # Mosaic and resample
    # far.data.country_biomass_compute(
    #     iso3=iso3,
    #     input_dir="data_raw",
    #     output_dir="data/emissions",
    #     proj=proj_aea)

    # # Remove GDAL temp directory
    # shutil.rmtree(temp_dir + "tmp_" + iso3)

    # Print country iso code
    print(iso3)


# Run country
for i in range(nctry):
    run_country(iso3[i])

# run_country(iso3[index_ctry])

# End
