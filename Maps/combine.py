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
import sys
import os
os.chdir("/home/gvieilledent/Code/forestatrisk-tropics/Maps/")
import subprocess
import forestatrisk as far
import matplotlib.pyplot as plt
from tif2cog import tif2cog

# Set PROJ_LIB
#PROJ_LIB = "/home/gvieilledent/.pyenv/versions/miniconda3-latest/envs/conda-far/share/proj"
#os.environ["PROJ_LIB"] = PROJ_LIB

# List of continents
continent = ["Africa", "America", "Asia"]
#ncont = len(continent)

# Loop on continent
#for i in range(ncont):

index_cont = int(sys.argv[1])-1

def run_combine(i):

    # Cont variable
    cont = continent[i]
    if cont == "America":
        cont_regex = "(America|Brazil)"
    else:
        cont_regex = cont

    # Make directory
    far.make_dir(os.path.expanduser("~/Code/forestatrisk-tropics/Maps/" + cont))
    os.chdir(os.path.expanduser("~/Code/forestatrisk-tropics/Maps/" + cont))

    # =======================
    # Combine country borders
    # =======================
    cmd = "find ~/nas/ -regextype posix-egrep -regex '.*" + cont_regex + ".*ctry_PROJ.shp$' \
    -exec ogr2ogr -update -append borders.shp {} \;"
    subprocess.call(cmd, shell=True)

    # # ===================
    # # Spatial probability
    # # ===================
    # # List of tif files
    # cmd = "find ~/nas/ -regextype posix-egrep -regex '.*" + cont_regex + ".*prob.tif$' > list_prob.txt"
    # subprocess.call(cmd, shell=True)
    # # COG
    # tif2cog(input_file_list="list_prob.txt", output_file="prob.tif", num_threads=8)
    # # Resample at 500m
    # subprocess.call("gdalwarp -overwrite -multi -wo 'NUM_THREADS=8' -wm 4096 -tap -r near -tr 500 500 -co 'COMPRESS=DEFLATE' \
    # -co 'PREDICTOR=2' -co 'BIGTIFF=YES' prob.tif prob_500m.tif", shell=True)
    # # Plot
    # prob = far.plot.prob("prob_500m.tif", output_file="prob.png", addo=False,
    #                      borders="borders.shp",
    #                      zoom=None, dpi=300,
    #                      lw=0.5, c="grey")
    # plt.close(prob)

    # ====================
    # Forest cover in 2050
    # ====================
    # List of tif files
    cmd = "find ~/nas/ -regextype posix-egrep -regex '.*" + cont_regex + ".*fcc_2050.tif$' > list_fcc_2050.txt"
    subprocess.call(cmd, shell=True)
    # COG
    tif2cog(input_file_list="list_fcc_2050.txt", output_file="fcc_2050.tif", num_threads=8)
    # Resample at 500m
    subprocess.call("gdalwarp -overwrite -multi -wo 'NUM_THREADS=8' -wm 4096 -tap -r near -tr 500 500 -co 'COMPRESS=LZW' \
    -co 'PREDICTOR=2' -co 'BIGTIFF=YES' fcc_2050.tif fcc_2050_500m.tif", shell=True)
    # # Plot
    # fcc = far.plot.fcc("fcc_2050_500m.tif", output_file="fcc_2050.png", addo=False,
    #                    borders="borders.shp",
    #                    zoom=None, dpi=300,
    #                    col_for=(34, 139, 34, 255),  # forest green
    #                    col_defor=(227, 26, 28, 255),  # red
    #                    lw=0.5, c="grey")
    # plt.close(fcc)

    # ==================================
    # Forest cover change 2000-2010-2019
    # ==================================
    cmd = "find ~/nas/ -regextype posix-egrep -regex '.*" + cont_regex + ".*data/fcc23.tif$' > list_fcc23.txt"
    subprocess.call(cmd, shell=True)
    # COG
    tif2cog(input_file_list="list_fcc23.txt", output_file="fcc23.tif", num_threads=8)
    # Resample at 500m
    subprocess.call("gdalwarp -overwrite -multi -wo 'NUM_THREADS=8' -wm 4096 -tap -r near -tr 500 500 -co 'COMPRESS=LZW' \
    -co 'PREDICTOR=2' -co 'BIGTIFF=YES' fcc23.tif fcc23_500m.tif", shell=True)
    # # Plot
    # fcc = far.plot.fcc("fcc23_500m.tif", output_file="fcc23.png", addo=False,
    #                    borders="borders.shp",
    #                    zoom=None, dpi=300,
    #                    col_for=(34, 139, 34, 255),  # forest green
    #                    col_defor=(227, 26, 28, 255),  # red
    #                    lw=0.5, c="grey")
    # plt.close(fcc)

#
run_combine(index_cont)

# End
