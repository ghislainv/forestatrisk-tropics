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
# 2. Google Cloud Storage: gcloud init

# Imports
import os
import forestatrisk as far
import pandas as pd
import pkg_resources
import psutil
import multiprocessing as mp
import matplotlib.pyplot as plt
os.chdir("/home/gvieilledent/Code/forestatrisk-tropics/Africa")
from run_plot import run_plot

# Set PROJ_LIB
PROJ_LIB = "/home/gvieilledent/.pyenv/versions/miniconda3-latest/envs/conda-far/share/proj"
os.environ["PROJ_LIB"] = PROJ_LIB

# Set working directory to nas
owd = "/share/nas2-amap/gvieilledent/Brazil"
os.chdir(owd)

# Country isocode
file_ctry_run = pkg_resources.resource_filename("forestatrisk",
                                                "data/ctry_run.csv")
data_ctry_run = pd.read_csv(file_ctry_run, sep=";", header=0)
iso3 = list(data_ctry_run.iso3[data_ctry_run.cont_run == "Brazil"])
iso3.sort()
# print(iso3)
# iso3 = ['AGO', 'BDI', 'BEN', 'CAF', 'CIV', 'CMR', 'COD', 'COG', 'COM',
#         'ETH', 'GAB', 'GHA', 'GIN', 'GMB', 'GNB', 'GNQ', 'KEN', 'LBR',
#         'MDG', 'MUS', 'MWI', 'MYT', 'NGA', 'REU', 'RWA', 'SEN',
#         'SLE', 'SSD', 'STP', 'TGO', 'TZA', 'UGA', 'ZMB']


# Number of cpu
total_cpu = psutil.cpu_count()
num_cpu = int(total_cpu * 0.75) if total_cpu > 2 else 1
# num_cpu = 4

# Function for multiprocessing
def run_country_plot(iso3):

    # Make new directory for country
    os.chdir(owd)
    far.make_dir(iso3)
    os.chdir(os.path.join(owd, iso3))

    # Model and Forecast
    run_plot()

    # Return country iso code
    return(iso3)

# For loop
#for i in iso3:
#    run_country(i)

# Parallel computation
pool = mp.Pool(processes=num_cpu)
results = [pool.apply_async(run_country_plot, args=(x,)) for x in iso3]
pool.close()
pool.join()
output = [p.get() for p in results]
print(output)

# End
