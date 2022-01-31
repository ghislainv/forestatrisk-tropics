#!/usr/bin/python
# -*- coding: utf-8 -*-

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ghislainv.github.io
# python_version  :>=3.0
# license         :GPLv3
# ==============================================================================

# Imports
import os
import subprocess
import sys

import forestatrisk as far

os.chdir(os.path.expanduser("~/Code/forestatrisk-tropics/Maps/"))
from leaflet_cog import leaflet_cog

# Set PROJ_LIB
# PROJ_LIB = "/home/ghislain/.pyenv/versions/miniconda3-latest/envs/conda-far/share/proj"
# os.environ["PROJ_LIB"] = PROJ_LIB

# List of continents
continent = ["Africa", "America", "Asia"]
ncont = len(continent)

# Directories
dataset = "jrc2020"
# rdir = "/home/forestatrisk-tropics/" + dataset
rdir = "/share/nas2-amap/gvieilledent/" + dataset

index_cont = int(sys.argv[1])-1


# run_combine
def run_combine(index_cont):

    # Cont variable
    cont = continent[index_cont]
    if cont == "America":
        cont_regex = "(America|Brazil)"
    else:
        cont_regex = cont
    # Make directory
    far.make_dir(rdir + "/Maps/" + cont)
    os.chdir(rdir + "/Maps/" + cont)

    # ===================
    # fcc123
    # ===================
    # List of tif files
    cmd = "find " + rdir + " -regextype posix-egrep -regex '.*" + \
        cont_regex + ".*fcc123.tif$' > list_tif.txt"
    subprocess.call(cmd, shell=True)
    # COG
    leaflet_cog(input_file_list="list_tif.txt",
                output_file="fcc123_epsg3857.tif", num_threads="ALL_CPUS")

    # ===================
    # Spatial probability
    # ===================
    # List of tif files
    cmd = "find " + rdir + " -regextype posix-egrep -regex '.*" + \
        cont_regex + ".*prob.tif$' > list_tif.txt"
    subprocess.call(cmd, shell=True)
    # COG
    leaflet_cog(input_file_list="list_tif.txt",
                output_file="prob_epsg3857.tif", num_threads="ALL_CPUS")

    # =========================================
    # Forest cover in 20XX - mean deforestation
    # =========================================

    # Dates
    dates_fut = [2050, 2100]
    ndates_fut = len(dates_fut)

    # Loop on dates
    for j in range(ndates_fut):
        # Date
        d = str(dates_fut[j])
        #
        if not os.path.isfile("fcc_" + d + "_epsg3857.tif"):
            # List of tif files
            cmd = "find " + rdir + " -regextype posix-egrep -regex '.*" + \
                cont_regex + ".*mean/fcc_" + d + ".tif$' > list_tif.txt"
            subprocess.call(cmd, shell=True)
            # COG
            leaflet_cog(input_file_list="list_tif.txt",
                        output_file="fcc_" + d + "_epsg3857.tif",
                        num_threads="ALL_CPUS")


# Run funtion
# for i in range(ncont):
#     run_combine(i)

run_combine(index_cont)

# EOF
