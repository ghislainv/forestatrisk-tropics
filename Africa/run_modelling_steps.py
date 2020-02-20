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


# run_modelling_steps
def run_modelling_steps(fcc_source="jrc"):

    # Make output directory
    far.make_dir("output_jrc")

    # ========================================================
    # Sample points
    # ========================================================

    dataset = far.sample(nsamp=10000, adapt=True,
                         Seed=1234, csize=10,
                         var_dir="data",
                         input_forest_raster="fcc23.tif",
                         output_file="output_jrc/sample.txt",
                         blk_rows=0)

    # To import data as pandas DataFrame if necessary
    # import pandas as pd
    # dataset = pd.read_table("output_jrc/sample.txt", delimiter=",")
    # dataset.head(5)

    # Descriptive statistics
    # Model formulas
    formula_1 = "fcc23 ~ dist_road + dist_town + dist_river + \
    dist_defor + dist_edge + altitude + slope - 1"
    # Standardized variables (mean=0, std=1)
    formula_2 = "fcc23 ~ scale(dist_road) + scale(dist_town) + \
    scale(dist_river) + scale(dist_defor) + scale(dist_edge) + \
    scale(altitude) + scale(slope) - 1"
    formulas = (formula_1, formula_2)

    # Remove NA from data-set (otherwise scale() and
    # model_binomial_iCAR doesn't work)
    dataset = dataset.dropna(axis=0)

    # Loop on formulas
    for f in range(len(formulas)):
        # Output file
        of = "output_jrc/correlation_" + str(f) + ".pdf"
        # Data
        y, data = dmatrices(formulas[f], data=dataset,
                            return_type="dataframe")
        # Plots
        figs = far.plot.correlation(y=y, data=data,
                                    plots_per_page=3,
                                    figsize=(7, 8),
                                    dpi=300,
                                    output_file=of)
        plt.close("all")

    # ========================================================
    # hSDM model
    # ========================================================

    # Set number of trials to one
    dataset["trial"] = 1

    # Spatial cells for spatial-autocorrelation
    nneigh, adj = far.cellneigh(raster="data/fcc23.tif", csize=10, rank=1)

    # List of variables
    variables = ["C(pa)", "scale(altitude)", "scale(slope)",
                 "scale(dist_defor)", "scale(dist_edge)", "scale(dist_road)",
                 "scale(dist_town)", "scale(dist_river)"]
    variables = np.array(variables)

    # Run model while there is non-significant variables
    var_remove = True
    while(np.any(var_remove)):

        # Formula
        right_part = " + ".join(variables) + " + cell"
        left_part = "I(1-fcc23) + trial ~ "
        formula = left_part + right_part

        # Model
        mod_binomial_iCAR = far.model_binomial_iCAR(
            # Observations
            suitability_formula=formula, data=dataset,
            # Spatial structure
            n_neighbors=nneigh, neighbors=adj,
            # Chains
            burnin=1000, mcmc=1000, thin=1,
            # Starting values
            beta_start=-99)

        # Ecological and statistical significance
        effects = mod_binomial_iCAR.betas[1:]
        # MCMC = mod_binomial_iCAR.mcmc
        # CI_low = np.percentile(MCMC, 2.5, axis=0)[1:-2]
        # CI_high = np.percentile(MCMC, 97.5, axis=0)[1:-2]
        positive_effects = (effects >= 0)
        # zero_in_CI = ((CI_low * CI_high) <= 0)

        # Keeping only significant variables
        var_remove = positive_effects
        # var_remove = np.logical_or(positive_effects, zero_in_CI)
        var_keep = np.logical_not(var_remove)
        variables = variables[var_keep]

    # Re-run the model with longer MCMC and estimated initial values
    mod_binomial_iCAR = far.model_binomial_iCAR(
        # Observations
        suitability_formula=formula, data=dataset,
        # Spatial structure
        n_neighbors=nneigh, neighbors=adj,
        # Chains
        burnin=5000, mcmc=5000, thin=5,
        # Starting values
        beta_start=mod_binomial_iCAR.betas)

    # Summary
    print(mod_binomial_iCAR)
    # Write summary in file
    f = open("output_jrc/summary_hSDM.txt", "w")
    f.write(str(mod_binomial_iCAR))
    f.close()

    # Plot
    figs = mod_binomial_iCAR.plot(output_file="output_jrc/mcmc.pdf",
                                  plots_per_page=3,
                                  figsize=(9, 6),
                                  dpi=300)
    plt.close("all")

    # ========================================================
    # Interpolating spatial random effects
    # ========================================================

    # Spatial random effects
    rho = mod_binomial_iCAR.rho

    # Interpolate
    far.interpolate_rho(rho=rho, input_raster="data/fcc23.tif",
                        output_file="output_jrc/rho.tif",
                        csize_orig=10, csize_new=1)

    # ========================================================
    # Predicting spatial probability of deforestation
    # ========================================================

    # Update dist_edge and dist_defor between t2 and t3
    os.rename("data/dist_edge.tif", "data/dist_edge.tif.bak")
    os.rename("data/dist_defor.tif", "data/dist_defor.tif.bak")
    copy2("data/proj/dist_edge_proj.tif", "data/dist_edge.tif")
    copy2("data/proj/dist_defor_proj.tif", "data/dist_defor.tif")

    # Compute predictions
    far.predict_raster_binomial_iCAR(
        mod_binomial_iCAR, var_dir="data",
        input_cell_raster="output_jrc/rho.tif",
        input_forest_raster="data/forest/forest_t3.tif",
        output_file="output_jrc/prob.tif",
        blk_rows=128
    )

    # ========================================================
    # Mean annual deforestation rate (ha.yr-1)
    # ========================================================

    # Forest cover
    fc = list()
    for i in range(3):
        rast = "data/forest/forest_t" + str(i+1) + ".tif"
        val = far.countpix(input_raster=rast,
                           value=1)
        fc.append(val["area"])
    # Save results to disk
    f = open("output_jrc/forest_cover.txt", "w")
    for i in fc:
        f.write(str(i) + "\n")
    f.close()

    # Annual deforestation
    T = 9.0 if (fcc_source == "jrc") else 9.0
    annual_defor = (fc[1] - fc[2]) / T
    # Amount of deforestation (ha)
    defor_2030 = np.rint(annual_defor * 11)
    defor_2050 = np.rint(annual_defor * 31)

    # ========================================================
    # Predicting forest cover change
    # ========================================================

    # Compute future forest cover
    stats = far.deforest(input_raster="output_jrc/prob.tif",
                         hectares=defor_2050,
                         output_file="output_jrc/fcc_2050.tif",
                         blk_rows=128)

    # Save stats to disk with pickle
    pickle.dump(stats, open("output_jrc/stats.pickle", "wb"))

    # Plot histograms of probabilities
    fig_freq = far.plot.freq_prob(stats,
                                  output_file="output_jrc/freq_prob.png")
    plt.close(fig_freq)

    # Forest cover change with half deforestation
    stats = far.deforest(input_raster="output_jrc/prob.tif",
                         hectares=np.rint(defor_2050 / 2.0),
                         output_file="output_jrc/fcc_2050_half.tif",
                         blk_rows=128)

    # Forest cover change after 10 years
    stats = far.deforest(input_raster="output_jrc/prob.tif",
                         hectares=defor_2030,
                         output_file="output_jrc/fcc_2030.tif",
                         blk_rows=128)

    # ========================================================
    # Figures
    # ========================================================

    # Forest in 2019
    fig_forest = far.plot.forest("data/forest/forest_t3.tif",
                                 borders="data/ctry_PROJ.shp",
                                 output_file="output_jrc/forest_t3.png")
    plt.close(fig_forest)

    # Forest-cover change 2010-2019
    fig_fcc = far.plot.fcc("data/forest/fcc13.tif",
                           borders="data/ctry_PROJ.shp",
                           output_file="output_jrc/fcc13.png")
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
                             borders="data/ctry_PROJ.shp",
                             output_file="output_jrc/prob.png")
    plt.close(fig_prob)

    # Forest-cover change 2019-2050
    fig_fcc_2050 = far.plot.fcc("output_jrc/fcc_2050.tif",
                                borders="data/ctry_PROJ.shp",
                                output_file="output_jrc/fcc_2050.png")
    plt.close(fig_fcc_2050)

    # Forest-cover change 2019-2030
    fig_fcc_2030 = far.plot.fcc("output_jrc/fcc_2030.tif",
                                borders="data/ctry_PROJ.shp",
                                output_file="output_jrc/fcc_2030.png")
    plt.close(fig_fcc_2030)

    # Forest-cover change 2019-2050 with half deforestation
    fig_fcc_2050_half = far.plot.fcc("output_jrc/fcc_2050_half.tif",
                                     borders="data/ctry_PROJ.shp",
                                     output_file="output_jrc/fcc_2050_half.png")
    plt.close(fig_fcc_2050_half)

# End
