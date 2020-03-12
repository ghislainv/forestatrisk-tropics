#!/usr/bin/python

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ghislainv.github.io
# python_version  :2.7
# license         :GPLv3
# ==============================================================================

# Imports
import os
import forestatrisk as far
import matplotlib.pyplot as plt


# run_modelling_steps
def run_plot():

    if os.path.isfile("output_jrc/prob.tif"):
        os.remove("output_jrc/prob.png")
        os.remove("output_jrc/prob.tif.ovr")
        fig_prob = far.plot.prob("output_jrc/prob.tif", addo=True,
                                 borders="data/ctry_PROJ.shp",
                                 output_file="output_jrc/prob.png")
        plt.close(fig_prob)

# End
