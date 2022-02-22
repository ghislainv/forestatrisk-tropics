#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ghislainv.github.io
# python_version  :2.7
# license         :GPLv3
# ==============================================================================

# Import
import os
import pkg_resources

import forestatrisk as far
import numpy as np
import pandas as pd

# Original working directory
owd = "/share/nas2-amap/gvieilledent/jrc2020/Brazil"
# owd = "/home/forestatrisk-tropics/jrc2020/Brazil"

# Country isocode for Brazil
file_ctry_run = pkg_resources.resource_filename("forestatrisk",
                                                "data/ctry_run.csv")
data_ctry_run = pd.read_csv(file_ctry_run, sep=";", header=0)
iso3 = data_ctry_run.loc[data_ctry_run["cont_run"] == "Brazil",
                         "iso3"].tolist()
# Sort iso3 (checking)
iso3.sort()
nctry = len(iso3)  # 26

# Loop on states to estimate the forest cover
for i in range(nctry):
    # Message
    print(iso3[i])
    # Directory
    os.chdir(os.path.join(owd, iso3[i]))
    far.make_dir("output")
    if os.path.exists("output/forest_cover.txt") is not True:
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
    fcc_BRA.iloc[i, 1:] = fc.iloc[:, 0].values.tolist()  # fc

# Dates for future predictions
dates_fut = ["2030", "2035", "2040", "2050", "2055", "2060",
             "2070", "2080", "2085", "2090", "2100", "2110",
             "2120", "2130", "2140", "2150"]
ndates_fut = len(dates_fut)

# Deforestation diffusion parameters
forest_t0 = np.array(fcc_BRA["for2020"])
t0 = 2020

# Deforestation estimates with uncertainty
f = os.path.expanduser("~/Code/forestatrisk-tropics/"
                       "Intensity/output/d_uncertainty.csv")
d_est = pd.read_csv(f)

# Take values for Brazil and **sort by iso3**
d_est_BRA = d_est[d_est["area_ctry"] == "Brazil"].sort_values("iso3")
# Check order
# d_est_BRA["iso3"].tolist() == iso3

# Scenarios
scenarios = ["mean", "min", "max"]
nscen = len(scenarios)

# Loop on scenarios
for k in range(nscen):

    scen = scenarios[k]

    # Copy fcc_BRA
    fcc_BRA_scen = fcc_BRA

    # Annual deforestation
    andef = d_est_BRA["d_" + scen].values

    # Loop on dates
    for i in range(ndates_fut):
        t = int(dates_fut[i])
        defor_diff = far.deforest_diffusion(
            forest_t0=forest_t0,
            t0=t0,
            annual_defor=andef.astype(float),
            t=t)
        fcc_BRA_scen["for" + dates_fut[i]] = defor_diff["forest_t"]
        fcc_BRA_scen["defor" + dates_fut[i]] = defor_diff["defor_t0_t"]

    # No more forest
    defor_diff_t_nofor = far.deforest_diffusion_t_nofor(
        forest_t0=forest_t0,
        t0=t0,
        annual_defor=andef.astype(float))
    fcc_BRA_scen["ny_zerof"] = defor_diff_t_nofor["ny"]
    fcc_BRA_scen["yr_zerof"] = defor_diff_t_nofor["year"]

    # Save results
    os.chdir(owd)
    ofile = "fcc_BRA_jrc_{}.csv".format(scen)
    fcc_BRA_scen.to_csv(ofile, header=True, index=False)

# EOF
