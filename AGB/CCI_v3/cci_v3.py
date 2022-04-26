#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ecology.ghislainv.fr
# python_version  :>=3
# license         :GPLv3
# ==============================================================================

# Standard library imports
import os
from ftplib import FTP
from glob import glob  # To explore files in a folder

# Third party imports
from osgeo import gdal


# ftp_download
def ftp_download(x):
    print("Downloading =>", x)
    fhandle = open(x, 'wb')
    ftp.retrbinary('RETR ' + x, fhandle.write)
    fhandle.close()


# progress bar
def progress_callback(complete, message, callback_data):
    '''Emit progress report in numbers for 10% intervals and dots for 3%'''
    if int(complete*100) % 10 == 0:
        print(f'{complete*100:.0f}', end='', flush=True)
    elif int(complete*100) % 3 == 0:
        print(f'{callback_data}', end='', flush=True)


# FTP
ftp = FTP("anon-ftp.ceda.ac.uk")
ftp.login()
ftp.cwd("/neodc/esacci/biomass/data/agb/maps/v3.0/geotiff/2010")
ftp_flist = ftp.nlst()

# File names
post_url = "_ESACCI-BIOMASS-L4-AGB-MERGED-100m-2010-fv3.0.tif"
lat_list = ["S30", "S20", "S10", "S00", "N00", "N10", "N20", "N30", "N40"]
lon_list_w = ["W{:03d}".format(x) for x in range(180, 0, -10)]
lon_list_e = ["E{:03d}".format(x) for x in range(0, 180, 10)]
lon_list = lon_list_w + lon_list_e

# Download
for lat in lat_list:
    for lon in lon_list:
        f = lat + lon
        fname = f + post_url
        if fname in ftp_flist and not os.path.isfile(fname):
            ftp_download(fname)

# Build VRT
file_list = glob("*_ESACCI-BIOMASS-L4-AGB-MERGED-100m-2010-fv3.0.tif")
file_vrt = "biomass_cci.vrt"  # path to vrt to build
gdal.BuildVRT(file_vrt, file_list)

# Build COG
ofile = "biomass_cci.tif"
param = gdal.TranslateOptions(
    options=["overwrite"],
    format="COG",
    creationOptions=["COMPRESS=LZW",
                     "PREDICTOR=YES",
                     "BIGTIFF=YES"],
    callback=progress_callback,
    callback_data=".")
gdal.Translate(ofile, file_vrt, options=param)

# Data cleaning
os.remove(file_vrt)
for i in file_list:
    os.remove(i)

# EOF
