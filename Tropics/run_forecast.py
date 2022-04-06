#!/usr/bin/python

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ghislainv.github.io
# python_version  :2.7
# license         :GPLv3
# ==============================================================================

import os
import re

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import forestatrisk as far


# run_forecast
def run_forecast(iso3, fcc_source="jrc"):
    """Forecast the forest cover in the future.

    Forecast the forest cover in the future following a
    business-as-usual scenario and estimate associated carbone
    emissions.

    :param iso3: Country ISO 3166-1 alpha-3 code. This is used to
    handle exceptions in the modelling process for specific countries.

    :param fcc_source: Source of forest cover data. Either "jrc" or
    "gfc".

    """

    # Make output directory
    far.make_dir("output")

    # ========================================================
    # Deforestation scenarios
    # ========================================================

    # Check if Brazil
    p = re.compile("BRA-.*")
    m = p.match(iso3)

    # Scenarios
    scenarios = ["mean", "min", "max"]
    nscen = len(scenarios)

    # Deforestation estimates with uncertainty
    f = os.path.expanduser("~/Code/forestatrisk-tropics/"
                           "Intensity/output/d_uncertainty.csv")
    d_est = pd.read_csv(f)

    # Loop on scenarios
    for k in range(nscen):

        # ---------------------------------------------------------
        # Preparing
        # --------------------------------------------------------

        scen = scenarios[k]

        # Make directory
        far.make_dir("output/" + scen)

        # Annual deforested area (ha.yr-1)
        annual_defor = d_est.loc[d_est["iso3"] == iso3,
                                 "d_" + scen].values[0]

        # For Brazil
        if m is not None:
            if (fcc_source == "jrc"):
                ifile = "../fcc_BRA_jrc_{}.csv".format(scen)
                fcc_BRA = pd.read_csv(ifile)
            else:
                ifile = "../fcc_BRA_gfc_{}.csv".format(scen)
                fcc_BRA = pd.read_csv(ifile)

        # Dates and time intervals
        dates_fut = ["2030", "2035", "2040", "2050", "2055", "2060",
                     "2070", "2080", "2085", "2090", "2100", "2110"]
        ndates_fut = len(dates_fut)
        ti = [10, 15, 20, 30, 35, 40, 50, 60, 65, 70, 80, 90]

        # --------------------------------------------------------
        # Predicting forest cover change
        # --------------------------------------------------------

        # # Loop on dates
        # for i in range(ndates_fut):
        #     # Amount of deforestation (ha)
        #     if m is not None:  # For Brazil
        #         defor = fcc_BRA.loc[fcc_BRA["iso3"] == iso3,
        #                             "defor" + dates_fut[i]].values[0]
        #     else:
        #         defor = np.rint(annual_defor * ti[i])
        #     # Compute future forest cover
        #     ofile = "output/{}/fcc_{}.tif".format(scen, dates_fut[i])
        #     far.deforest(
        #         input_raster="output/prob.tif",
        #         hectares=defor,
        #         output_file=ofile,
        #         blk_rows=128)

        # --------------------------------------------------------
        # Carbon emissions
        # --------------------------------------------------------

        # Create dataframe
        dpast = ["2020"]
        dpast.extend(dates_fut)
        C_df = pd.DataFrame(
            {"date": dpast, "C": np.repeat(-99, ndates_fut + 1)},
            columns=["date", "C"])
        # Loop on date
        for i in range(ndates_fut):
            ifile = "output/{}/fcc_{}.tif".format(scen, dates_fut[i])
            carbon = far.emissions(
                input_stocks="data/emissions/biomass_whrc.tif",
                input_forest=ifile)
            C_df.loc[C_df["date"] == dates_fut[i], ["C"]] = carbon
        # Past emissions
        carbon = far.emissions(input_stocks="data/emissions/biomass_whrc.tif",
                               input_forest="data/fcc23.tif")
        C_df.loc[C_df["date"] == dpast[0], ["C"]] = carbon
        # Save dataframe
        ofile = "output/{}/C_emissions_whrc.csv".format(scen)
        C_df.to_csv(ofile, header=True, index=False)

        # --------------------------------------------------------
        # Figures
        # --------------------------------------------------------

        # # Projected future forest-cover change
        # for i in range(ndates_fut):
        #     ifile = "output/{}/fcc_{}.tif".format(scen, dates_fut[i])
        #     ofile = "output/{}/fcc_{}.png".format(scen, dates_fut[i])
        #     fig_fcc = far.plot.fcc(
        #         ifile,
        #         maxpixels=1e8,
        #         borders="data/ctry_PROJ.shp",
        #         output_file=ofile)
        #     plt.close(fig_fcc)

# End
