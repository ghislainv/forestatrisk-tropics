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
import sys
import os
# os.chdir("/home/ghislain/Code/forestatrisk-tropics/Maps/")
os.chdir("/home/gvieilledent/Code/forestatrisk-tropics/Maps/")
import subprocess
import forestatrisk as far
import matplotlib.pyplot as plt
from tif2cog import tif2cog

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

    # =======================
    # Combine roads
    # =======================
    # Combine
    cmd = "find " + rdir + " -regextype posix-egrep -regex '.*" + cont_regex + ".*/data/roads_PROJ.shp$' \
    -exec ogr2ogr -update -append -nlt MULTILINESTRING roads.gpkg {} \;"
    subprocess.call(cmd, shell=True)
    # Simplify and reproject
    cmd = "ogr2ogr -overwrite -t_srs EPSG:3395 -nlt MULTILINESTRING roads_simp.gpkg roads.gpkg -simplify 1000"
    subprocess.call(cmd, shell=True)

    # =======================
    # Combine protected areas
    # =======================
    # Combine
    cmd = "find " + rdir + " -regextype posix-egrep -regex '.*" + cont_regex + ".*/data/pa_PROJ.shp$' \
    -exec ogr2ogr -update -append -nlt MULTIPOLYGON pa.gpkg {} \;"
    subprocess.call(cmd, shell=True)
    # Simplify
    cmd = "ogr2ogr -overwrite -t_srs EPSG:3395 -nlt MULTIPOLYGON pa_simp.gpkg pa.gpkg -simplify 1000"
    subprocess.call(cmd, shell=True)

    # =======================
    # Combine country borders
    # =======================
    # Combine
    cmd = "find " + rdir + " -regextype posix-egrep -regex '.*" + cont_regex + ".*/data/ctry_PROJ.shp$' \
    -exec ogr2ogr -update -append -nlt MULTIPOLYGON borders.gpkg {} \;"
    subprocess.call(cmd, shell=True)
    # Simplify
    cmd = "ogr2ogr -overwrite -t_srs EPSG:3395 -nlt MULTIPOLYGON borders_simp.gpkg borders.gpkg -simplify 1000"
    subprocess.call(cmd, shell=True)

    # ===================
    # fcc123
    # ===================
    # List of tif files
    cmd = "find " + rdir + " -regextype posix-egrep -regex '.*" + cont_regex + ".*fcc123.tif$' > list_tif.txt"
    subprocess.call(cmd, shell=True)
    # COG
    tif2cog(input_file_list="list_tif.txt", output_file="fcc123.tif", num_threads=8)
    # Resample at 500m
    subprocess.call("gdalwarp -overwrite -multi -wo 'NUM_THREADS=8' -wm 4096 -tap -r near -tr 500 500 -co 'COMPRESS=DEFLATE' \
    -co 'PREDICTOR=2' -co 'BIGTIFF=YES' fcc123.tif fcc123_500m.tif", shell=True)
    # Plot
    fcc123 = far.plot.fcc123("fcc123_500m.tif", output_file="fcc123.png",
                             maxpixels=1e13,
                             borders="borders_simp.gpkg",
                             zoom=None, dpi=300,
                             lw=0.5, c="grey")
    plt.close(fcc123)

    # ===================
    # Spatial probability
    # ===================
    # List of tif files
    cmd = "find " + rdir + " -regextype posix-egrep -regex '.*" + cont_regex + ".*prob.tif$' > list_tif.txt"
    subprocess.call(cmd, shell=True)
    # COG
    tif2cog(input_file_list="list_tif.txt", output_file="prob.tif", num_threads=8)
    # Resample at 500m
    subprocess.call("gdalwarp -overwrite -multi -wo 'NUM_THREADS=8' -wm 4096 -t_srs EPSG:3395 \
    -tap -r near -tr 500 500 -co 'COMPRESS=DEFLATE' \
    -co 'PREDICTOR=2' -co 'BIGTIFF=YES' prob.tif prob_500m.tif", shell=True)
    # Plot
    prob = far.plot.prob("prob_500m.tif", output_file="prob.png",
                         maxpixels=1e13,
                         borders="borders_simp.gpkg",
                         zoom=None, dpi=300,
                         lw=0.5, c="grey")
    plt.close(prob)

    # ====================
    # Forest cover in 20XX
    # ====================

    # Dates
    dates_fut = [2030, 2035, 2040, 2050, 2055, 2060,
                 2070, 2080, 2085, 2090, 2100]
    ndates_fut = len(dates_fut)

    # Loop on dates
    for j in range(ndates_fut):
        # Date
        d = str(dates_fut[j])
        #
        if not os.path.isfile("fcc_" + d + ".tif"):
            # List of tif files
            cmd = "find " + rdir + " -regextype posix-egrep -regex '.*" + cont_regex + ".*fcc_" + d + ".tif$' > list_tif.txt"
            subprocess.call(cmd, shell=True)
            # COG
            tif2cog(input_file_list="list_tif.txt", output_file="fcc_" + d + ".tif", num_threads=8)
            # Resample at 500m
            subprocess.call("gdalwarp -overwrite -multi -wo 'NUM_THREADS=8' -wm 4096 -t_srs EPSG:3395 \
            -tap -r near -tr 500 500 -co 'COMPRESS=LZW' \
            -co 'PREDICTOR=2' -co 'BIGTIFF=YES' fcc_" + d + ".tif fcc_" + d + "_500m.tif", shell=True)
            # Plot
            fcc = far.plot.fcc("fcc_" + d + "_500m.tif", output_file="fcc_" + d + ".png",
                               maxpixels=1e13,
                               borders="borders_simp.gpkg",
                               zoom=None, dpi=300,
                               lw=0.5, c="grey")
            plt.close(fcc)

# Run funtion
# Loop on continent
#for index_cont in range(ncont):
run_combine(index_cont)


# # ===================
# # Resample for Maxime
# # ===================

# subprocess.call("gdalinfo AnthropogenicPressure.tif", shell=True)
# subprocess.call("gdalinfo prob.tif", shell=True)
# subprocess.call("gdalwarp -overwrite -multi -wo 'NUM_THREADS=8' -wm 4096 -r average \
# -ot UInt16 -t_srs EPSG:4326 \
# -te 7.9999248 -6.0000216 29.9999160 7.4999730 -tr 0.00833333 0.00833333 \
# -co 'COMPRESS=DEFLATE' -co 'PREDICTOR=2' -co 'BIGTIFF=YES' \
# prob.tif prob_10km.tif", shell=True)

# EOF
