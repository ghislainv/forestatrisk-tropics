#!/usr/bin/env python

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ghislainv.github.io
# python_version  :2.7
# license         :GPLv3
# ==============================================================================

# Import
import sys
import os
import subprocess
import pkg_resources
import pandas as pd
import numpy as np
import forestatrisk as far

# Original working directory
owd = "/share/nas2-amap/gvieilledent/jrc2020/Brazil"

# Country isocode for Brazil
file_ctry_run = pkg_resources.resource_filename("forestatrisk",
                                                "data/ctry_run.csv")
data_ctry_run = pd.read_csv(file_ctry_run, sep=";", header=0)
iso3 = data_ctry_run.loc[data_ctry_run["cont_run"]=="Brazil", "iso3"].tolist()
nctry = len(iso3)  # 26

# Loop on states to estimate the forest cover
for i in range(nctry):
    # Message
    print(iso3[i])
    # Directory
    os.chdir(os.path.join(owd, iso3[i]))
    far.make_dir("output")
    # Forest cover
    fc = list()
    dates = ["t1", "2005", "t2", "2015", "t3"]
    ndates = len(dates)
    for i in range(ndates):
        rast = "data/forest/forest_" + dates[i] + ".tif"
        val = far.countpix(input_raster=rast,
                           value=1)
        fc.append(val["area"])  # area in ha
    # Save results to disk
    f = open("output/forest_cover.txt", "w")
    for i in fc:
        f.write(str(i) + "\n")
    f.close()

# Create data-frame to store results
fcc_BRA = pd.DataFrame({"iso3": iso3},
                       columns=["iso3", "for2000", "for2005",
                                "for2010", "for2015", "for2020"])

# Fill in the data-frame
for i in range(nctry):
    # Directory
    os.chdir(os.path.join(owd, iso3[i]))
    # Read forest cover file
    fc = pd.read_csv("output/forest_cover.txt", header=None)
    # Fill in the table
    fcc_BRA.loc[i, 1:] = fc.loc[:, 0].values.tolist() # fc

# Dates for future predictions
dates_fut = ["2030", "2035", "2040", "2050", "2055", "2060", "2070", "2080", "2085", "2090", "2100"]
ndates_fut = len(dates_fut)

# Deforestation diffusion parameters
forest_t0 = np.array(fcc_BRA["for2020"])
T = 10.0
andef = np.array((fcc_BRA["for2010"] - fcc_BRA["for2020"]) / T)
t0 = 2020

# Loop on dates
for i in range(ndates_fut):  
    t = int(dates_fut[i])
    defor_diff = far.deforest_diffusion(forest_t0=forest_t0,
                                        t0=t0,
                                        annual_defor=andef,
                                        t=t)
    fcc_BRA["for" + dates_fut[i]] = defor_diff["forest_t"]
    fcc_BRA["defor" + dates_fut[i]] = defor_diff["defor_t0_t"]

# No more forest
defor_diff_t_nofor = far.deforest_diffusion_t_nofor(
    forest_t0=forest_t0,
    t0=t0,
    annual_defor=andef)
fcc_BRA["ny_zerof"] = defor_diff_t_nofor["ny"]
fcc_BRA["yr_zerof"] = defor_diff_t_nofor["year"]

# Save results
os.chdir(owd)
fcc_BRA.to_csv("fcc_BRA_jrc.csv", header=True, index=False)

# EOF
