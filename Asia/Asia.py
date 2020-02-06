#!/usr/bin/python
# -*- coding: utf-8 -*-

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ghislainv.github.io
# python_version  :2.7
# license         :GPLv3
# ==============================================================================

# PRIOR TO EXECUTING THE FOLLOWING SCRIPT, AUTHENTICATE TO
# 1. Google Earth Engine: earthengine authenticate
# 2. Google Drive with rclone: https://rclone.org/drive/

import os
import forestatrisk as far
import pandas as pd
import pkg_resources
import psutil
import multiprocessing as mp
from run_modelling_steps import run_modelling_steps

# Set PROJ_LIB
os.environ["PROJ_LIB"] = "/home/ghislain/miniconda3/envs/forestatrisk/share/proj"

# Set wd
#os.chdir("/home/ghislain/Code/forestatrisk-tropics/Asia")

# List of countries to process
# countries = ["Bangladesh", "Bhutan", "Cambodia", "Indonesia", "India",
#              "Lao People's Democratic Republic", "Malaysia", "Myanmar",
#              "Nepal", "New Caledonia", "Papua New Guinea", "Philippines",
#              "Sri Lanka", "Viet Nam", "Thailand"]

countries = ["Bangladesh", "Bhutan", "Cambodia", "India",
             "Lao People's Democratic Republic", "Malaysia", "Myanmar",
             "Nepal", "New Caledonia", "Papua New Guinea", "Philippines",
             "Sri Lanka", "Viet Nam", "Thailand"]
             
# Number of countries
nctry = len(countries)

# Data-frame of country codes
file_countrycode = pkg_resources.resource_filename("forestatrisk",
                                                   "data/countrycode.csv")
data_countrycode = pd.read_csv(file_countrycode, sep=";", header=0)

# Get iso3c from country name
iso3 = list()
for i in range(nctry):
    code = data_countrycode.iso3c[
        data_countrycode["country.name.en"] == countries[i]]
    iso3.append(code.iloc[0])
iso3.sort()

# Only some countries for test
#iso3 = ["BGD","BTN","IDN", "KHM","LAO", "LKA", "MMR", "MYS", "PNG", "THA", "VNM"]
#iso3 = ["VNM"]
nctry = len(iso3)

# Projection for Asia (World Mercator)
proj_asia = "EPSG:3395"

# Original working directory
owd = os.getcwd()

# Number of cpu
total_cpu = psutil.cpu_count()
num_cpu = int(total_cpu * 0.75) if total_cpu > 2 else 1
# num_cpu = 3


# Function for multiprocessing
def run_country(iso3):
    
    # Make new directory for country
    os.chdir(owd)
    far.make_dir(iso3)
    os.chdir(os.path.join(owd, iso3))
    
    # Data
    far.data.country(iso3=iso3, monthyear="Jan2020",
                     proj=proj_asia,
                     data_country=True,
                     data_forest=False,
                     keep_data_raw=True,
                     fcc_source="jrc",
                     gdrive_remote_rclone="gdrive_gv",
                     gdrive_folder="GEE-forestatrisk-tropics")
    
    # # Model and Forecast
    # run_modelling_steps(fcc_source="roadless")
    
    # Return country iso code
    return(iso3)

# For loop
for i in iso3:
    print("GEE for country: " + i + "\n")
    far.make_dir(i + "/data_raw")
    far.country_forest_gdrive(
        iso3=i, proj="EPSG:3395",
        output_dir=i + "/data_raw",
        keep_dir=True,
        fcc_source="jrc", perc=50,
        gdrive_remote_rclone="gdrive_gv",
        gdrive_folder="GEE-forestatrisk-tropics"
    )
   #run_country(i)

# # Parallel computation
# pool = mp.Pool(processes=num_cpu)
# results = [pool.apply_async(run_country, args=(x,)) for x in iso3]
# pool.close()
# pool.join()
# output = [p.get() for p in results]
# print(output)

# Combine results

# Combine country borders
# os.system("find -type f -name *ctry_PROJ.shp \
# -exec ogr2ogr -update -append borders.shp {} \;")

# Spatial probability
# os.system("find -type f -name *prob.tif > list_prob.txt")
# os.system("gdalbuildvrt -input_file_list list_prob.txt prob.vrt")
# os.system("gdal_translate -co 'COMPRESS=LZW' -co 'PREDICTOR=2' -co 'BIGTIFF=YES' \
# prob.vrt prob.tif")
# Build overview
# os.system("gdaladdo -ro -r nearest prob.tif 16")
# Plot
# dfp.plot.prob("prob.tif", output_file="prob.png",
#               borders="borders.shp", zoom=None, dpi=300,
#               lw=0.5, c="grey")

# Forest cover in 2050
# os.system("find -type f -name *fcc_40yr.tif > list_fcc_40yr.txt")
# os.system("gdalbuildvrt -input_file_list list_fcc_40yr.txt fcc_40yr.vrt")
# os.system("gdal_translate -co 'COMPRESS=LZW' -co 'PREDICTOR=2' -co 'BIGTIFF=YES' \
# fcc_40yr.vrt fcc_40yr.tif")
# Build overview
# os.system("gdaladdo -ro -r nearest --config COMPRESS_OVERVIEW LZW \
# --config PREDICTOR_OVERVIEW 2 \
# --config BIGTIFF_OVERVIEW YES \
# fcc_40yr.tif 16")
# Plot
# dfp.plot.fcc("fcc_40yr.tif", output_file="fcc_40yr.png",
#              borders="borders.shp", overview=False, zoom=None, dpi=300,
#              lw=0.5, c="grey")

# Upload results to Google Cloud Storage with gsutil
# os.system("gsutil -o GSUtil:parallel_composite_upload_threshold=150M \
# cp fcc_40yr.tif gs://deforestprob/output_roadless/fcc_40yr.tif")
# os.system("gsutil cp fcc_40yr.png \
# gs://deforestprob/output_roadless/fcc_40yr.png")
# os.system("gsutil -o GSUtil:parallel_composite_upload_threshold=150M \
# cp prob.tif gs://deforestprob/output_roadless/prob.tif")
# os.system("gsutil cp prob.png gs://deforestprob/output_roadless/prob.png")

# End
