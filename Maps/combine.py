#!/usr/bin/python
# -*- coding: utf-8 -*-

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ghislainv.github.io
# python_version  :2.7
# license         :GPLv3
# ==============================================================================

# Imports
import os
import subprocess
import forestatrisk as far
import matplotlib.pyplot as plt

# Set PROJ_LIB
#PROJ_LIB = "/home/gvieilledent/.pyenv/versions/miniconda3-latest/envs/conda-far/share/proj"
#os.environ["PROJ_LIB"] = PROJ_LIB

# List of continents
continent = ["Africa", "America", "Asia"]
ncont = len(continent)

# Working directory
owd = os.path.expanduser("~/Code/forestatrisk-tropics/Maps")
os.chdir(owd)

# Loop on continent
for i in range(ncont):

    # Cont variable
    cont = continent[i]
    if cont == "America":
        cont_regex = "(America|Brazil)"
    else:
        cont_regex = cont

    # Make directory
    far.make_dir(cont)
    os.chdir(os.path.expanduser("~/Code/forestatrisk-tropics/Maps/" + cont))

    # =======================
    # Combine country borders
    # =======================
    cmd = "find ~/nas/ -regextype posix-egrep -regex '.*" + cont_regex + ".*ctry_PROJ.shp$' \
    -exec ogr2ogr -update -append borders.shp {} \;"
    subprocess.call(cmd, shell=True)

    # ===================
    # Spatial probability
    # ===================
    cmd = "find ~/nas/ -regextype posix-egrep -regex '.*" + cont_regex + ".*prob.tif$' > list_prob.txt"
    subprocess.call(cmd, shell=True)
    subprocess.call("gdalbuildvrt -input_file_list list_prob.txt prob.vrt", shell=True)
    cmd = "gdal_translate -co 'COMPRESS=LZW' -co 'PREDICTOR=2' -co 'BIGTIFF=YES' \
    prob.vrt prob.tif"
    subprocess.call(cmd, shell=True)
    # Resample at 500m
    subprocess.call("gdalwarp -tap -r near -tr 500 500 -co 'COMPRESS=LZW' \
    -co 'PREDICTOR=2' -co 'BIGTIFF=YES' prob.tif prob_500m.tif", shell=True)
    # Plot
    prob = far.plot.prob("prob_500m.tif", output_file="prob.png", addo=False,
                borders="borders.shp", zoom=None, dpi=300,
                lw=0.5, c="grey")
    plt.close(prob)

    # ====================
    # Forest cover in 2050
    # ====================
    cmd = "find ~/nas/ -regextype posix-egrep -regex '.*" + cont_regex + ".*fcc_2050.tif$' > list_fcc_2050.txt"
    subprocess.call(cmd, shell=True)
    subprocess.call("gdalbuildvrt -input_file_list list_fcc_2050.txt fcc_2050.vrt", shell=True)
    cmd = "gdal_translate -co 'COMPRESS=LZW' -co 'PREDICTOR=2' -co 'BIGTIFF=YES' \
    fcc_2050.vrt fcc_2050.tif"
    subprocess.call(cmd, shell=True)
    # Resample at 500m
    subprocess.call("gdalwarp -tap -r near -tr 500 500 -co 'COMPRESS=LZW' \
    -co 'PREDICTOR=2' -co 'BIGTIFF=YES' fcc_2050.tif fcc_2050_500m.tif", shell=True)
    # Plot
    fcc = far.plot.fcc("fcc_2050_500m.tif", output_file="fcc_2050.png", addo=False,
                       borders="borders.shp", zoom=None, dpi=300,
                       lw=0.5, c="grey")
    plt.close(fcc)
    # Plot for colorblindless
    fcc = far.plot.fcc("fcc_2050_500m.tif", output_file="fcc_2050_colblind.png",
                       borders="borders.shp", zoom=None, dpi=300,
                       col_for=(30, 136, 229, 255),  # blue
                       #col_for=(0, 77, 64, 255),  # green
                       col_defor=(255, 193, 7, 255),  # yellow
                       lw=0.5, c="grey")
    plt.close(fcc)

    # =============================
    # Forest cover change 2019-2050
    # =============================
    cmd = "find ~/nas/ -regextype posix-egrep -regex '.*" + cont_regex + ".*fcc23.tif$' > list_fcc23.txt"
    subprocess.call(cmd, shell=True)
    subprocess.call("gdalbuildvrt -input_file_list list_fcc_2050.txt fcc_2050.vrt", shell=True)
    cmd = "gdal_translate -co 'COMPRESS=LZW' -co 'PREDICTOR=2' -co 'BIGTIFF=YES' \
    fcc_2050.vrt fcc_2050.tif"
    subprocess.call(cmd, shell=True)
    # Resample at 500m
    subprocess.call("gdalwarp -tap -r near -tr 500 500 -co 'COMPRESS=LZW' \
    -co 'PREDICTOR=2' -co 'BIGTIFF=YES' fcc_2050.tif fcc_2050_500m.tif", shell=True)
    # Plot
    fcc = far.plot.fcc("fcc_2050_500m.tif", output_file="fcc_2050.png", addo=False,
                       borders="borders.shp", zoom=None, dpi=300,
                       lw=0.5, c="grey")
    plt.close(fcc)
    # Plot for colorblindless
    fcc = far.plot.fcc("fcc_2050_500m.tif", output_file="fcc_2050_colblind.png",
                       borders="borders.shp", zoom=None, dpi=300,
                       col_for=(30, 136, 229, 255),  # blue
                       #col_for=(0, 77, 64, 255),  # green
                       col_defor=(255, 193, 7, 255),  # yellow
                       lw=0.5, c="grey")
    plt.close(fcc)
    

    # Return to owd
    os.chdir(owd)

# End
