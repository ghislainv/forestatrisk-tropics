#!/usr/bin/python

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ghislainv.github.io
# python_version  :2.7
# license         :GPLv3
# ==============================================================================

import os
from shutil import copy2  # To copy files
import numpy as np
from patsy import dmatrices
import forestatrisk as far
import matplotlib.pyplot as plt
import pickle
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss
import pandas as pd


# run_modelling_steps
def run_modelling_steps(fcc_source="jrc"):

    # Make output directory
    far.make_dir("output_jrc")

    # ========================================================
    # Mean annual deforestation rate (ha.yr-1)
    # ========================================================

    # Forest cover
    fc = list()
    for i in range(3):
        rast = "data/forest/forest_t" + str(i+1) + ".tif"
        val = far.countpix(input_raster=rast,
                           value=1)
        fc.append(val["area"])  # area in ha
    # Save results to disk
    f = open("output_jrc/forest_cover.txt", "w")
    for i in fc:
        f.write(str(i) + "\n")
    f.close()

    # Annual deforestation
    T = 9.0 if (fcc_source == "jrc") else 9.0
    annual_defor = (fc[1] - fc[2]) / T

    # Dates and time intervals
    date = ["2035", "2050", "2055", "2085", "2100"]
    ndate = len(date)
    ti = [16, 31, 36, 66, 81]
    
    # ========================================================
    # Predicting forest cover change
    # ========================================================

    # Loop on dates
    for i in range(ndate):
        # Amount of deforestation (ha)
        defor = np.rint(annual_defor * ti[i])
        # Compute future forest cover
        stats = far.deforest(input_raster="output_jrc/prob.tif",
                             hectares=defor,
                             output_file="output_jrc/fcc_" + date[i] + ".tif",
                             blk_rows=128)
        # Save some stats if date = 2050
        if date[i] == "2050":
            # Save stats to disk with pickle
            pickle.dump(stats, open("output_jrc/stats.pickle", "wb"))
            # Plot histograms of probabilities
            fig_freq = far.plot.freq_prob(stats,
                                          output_file="output_jrc/freq_prob.png")
            plt.close(fig_freq)

    # ========================================================
    # Figures
    # ========================================================

    # Forest in 2019
    fig_forest = far.plot.forest("data/forest/forest_t3.tif",
                                 maxpixels=5000000,
                                 borders="data/ctry_PROJ.shp",
                                 output_file="output_jrc/forest_t3.png")
    plt.close(fig_forest)

    # Forest-cover change 2000-2019
    fig_fcc = far.plot.fcc("data/forest/fcc13.tif",
                           maxpixels=5000000,
                           borders="data/ctry_PROJ.shp",
                           output_file="output_jrc/fcc13.png")
    plt.close(fig_fcc)

    # Forest-cover change 2010-2019
    fig_fcc = far.plot.fcc("data/fcc23.tif",
                           maxpixels=5000000,
                           borders="data/ctry_PROJ.shp",
                           output_file="output_jrc/fcc23.png")
    plt.close(fig_fcc)

    # Original spatial random effects
    fig_rho_orig = far.plot.rho("output_jrc/rho_orig.tif",
                                borders="data/ctry_PROJ.shp",
                                output_file="output_jrc/rho_orig.png")
    plt.close(fig_rho_orig)

    # Interpolated spatial random effects
    fig_rho = far.plot.rho("output_jrc/rho.tif",
                           borders="data/ctry_PROJ.shp",
                           output_file="output_jrc/rho.png")
    plt.close(fig_rho)

    # Spatial probability of deforestation
    fig_prob = far.plot.prob("output_jrc/prob.tif",
                             maxpixels=5000000,
                             borders="data/ctry_PROJ.shp",
                             output_file="output_jrc/prob.png")
    plt.close(fig_prob)

    # Projected future forest-cover change
    for i in range(ndate):
        fig_fcc = far.plot.fcc("output_jrc/fcc_" + date[i] + ".tif",
                               maxpixels=5000000,
                               borders="data/ctry_PROJ.shp",
                               output_file="output_jrc/fcc_" + date[i] + ".png")
        plt.close(fig_fcc)


# Temp
os.chdir("/home/gvieilledent/nas/Africa/MDG/")
run_modelling_steps()

# End
