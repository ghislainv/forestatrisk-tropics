#!/usr/bin/python

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ghislainv.github.io
# python_version  :2.7
# license         :GPLv3
# ==============================================================================

import os
import numpy as np
from patsy import dmatrices
import forestatrisk as far
import matplotlib.pyplot as plt
import pickle


# run_modelling_steps
def run_modelling_steps(fcc_source="roadless"):

    # Make output directory
    far.make_dir("output_roadless")

    # ========================================================
    # Sample points
    # ========================================================

    dataset = far.sample(nsamp=10000, Seed=1234, csize=10,
                         var_dir="data",
                         input_forest_raster="fcc23.tif",
                         output_file="output_roadless/sample.txt",
                         blk_rows=0)

    # To import data as pandas DataFrame if necessary
    # import pandas as pd
    # dataset = pd.read_table("output_roadless/sample.txt", delimiter=",")
    # dataset.head(5)

    # Descriptive statistics
    # Model formulas
    formula_1 = "fcc23 ~ dist_road + dist_town + dist_river + \
    dist_defor + dist_edge + altitude + slope + aspect - 1"
    # Standardized variables (mean=0, std=1)
    formula_2 = "fcc23 ~ scale(dist_road) + scale(dist_town) + \
    scale(dist_river) + scale(dist_defor) + scale(dist_edge) + \
    scale(altitude) + scale(slope) + scale(aspect) - 1"
    formulas = (formula_1, formula_2)

    # Remove NA from data-set (otherwise scale() and
    # model_binomial_iCAR doesn't work)
    dataset = dataset.dropna(axis=0)

    # Loop on formulas
    for f in range(len(formulas)):
        # Output file
        of = "output_roadless/correlation_" + str(f) + ".pdf"
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
    f = open("output_roadless/summary_hSDM.txt", "w")
    f.write(str(mod_binomial_iCAR))
    f.close()

    # Plot
    figs = mod_binomial_iCAR.plot(output_file="output_roadless/mcmc.pdf",
                                  plots_per_page=3,
                                  figsize=(9, 6),
                                  dpi=300)
    plt.close("all")

    # ========================================================
    # Resampling spatial random effects
    # ========================================================

    # Spatial random effects
    rho = mod_binomial_iCAR.rho

    # Resample
    far.resample_rho(rho=rho, input_raster="data/fcc23.tif",
                     output_file="output_roadless/rho.tif",
                     csize_orig=10, csize_new=1)

    # ========================================================
    # Predicting spatial probability of deforestation
    # ========================================================

    # We assume dist_edge and dist_defor don't change between t2
    # and t3 (deforestation ~1%). No need to recompute them.

    # Rename aspect.tif in data directory to avoid NA where slope=0
    os.rename("data/aspect.tif", "data/aspect.tif.bak")

    # Compute predictions
    far.predict(mod_binomial_iCAR, var_dir="data",
                input_cell_raster="output_roadless/rho.tif",
                input_forest_raster="data/forest/forest_t3.tif",
                output_file="output_roadless/prob.tif",
                blk_rows=128)

    # Rename aspect.tif.bak
    os.rename("data/aspect.tif.bak", "data/aspect.tif")

    # ========================================================
    # Mean annual deforestation rate (ha.yr-1)
    # ========================================================

    # Forest cover
    fc = list()
    for i in range(4):
        rast = "data/forest/forest_t" + str(i) + ".tif"
        val = far.countpix(input_raster=rast,
                           value=1)
        fc.append(val["area"])
    # Save results to disk
    f = open("output_roadless/forest_cover.txt", "w")
    for i in fc:
        f.write(str(i) + "\n")
    f.close()

    # Annual deforestation
    T = 10.0 if (fcc_source == "roadless") else 9.0
    annual_defor = (fc[1] - fc[3]) / T
    # Amount of deforestation (ha)
    defor_10yr = np.rint(annual_defor * 10)
    defor_35yr = np.rint(annual_defor * 35)

    # ========================================================
    # Predicting forest cover change
    # ========================================================

    # Compute future forest cover
    stats = far.deforest(input_raster="output_roadless/prob.tif",
                         hectares=defor_35yr,
                         output_file="output_roadless/fcc_35yr.tif",
                         blk_rows=128)

    # Save stats to disk with pickle
    pickle.dump(stats, open("output_roadless/stats.pickle", "wb"))

    # Plot histograms of probabilities
    fig_freq = far.plot.freq_prob(stats,
                                  output_file="output_roadless/freq_prob.png")
    plt.close(fig_freq)
    
    # Forest cover change with half deforestation
    stats = far.deforest(input_raster="output_roadless/prob.tif",
                         hectares=np.rint(defor_35yr / 2.0),
                         output_file="output_roadless/fcc_35yr_half.tif",
                         blk_rows=128)
                         
    # Forest cover change after 10 years
    stats = far.deforest(input_raster="output_roadless/prob.tif",
                         hectares=defor_10yr,
                         output_file="output_roadless/fcc_10yr.tif",
                         blk_rows=128)

    # ========================================================
    # Figures
    # ========================================================

    # Forest in 2015
    fig_forest = far.plot.forest("data/forest/forest_t3.tif",
                                 borders="data/ctry_PROJ.shp",
                                 output_file="output_roadless/forest_t3.png")
    plt.close(fig_forest)

    # Forest-cover change 2005-2015
    fig_fcc = far.plot.fcc("data/forest/fcc13.tif",
                           borders="data/ctry_PROJ.shp",
                           output_file="output_roadless/fcc13.png")
    plt.close(fig_fcc)

    # Original spatial random effects
    fig_rho_orig = far.plot.rho("output_roadless/rho_orig.tif",
                                borders="data/ctry_PROJ.shp",
                                output_file="output_roadless/rho_orig.png")
    plt.close(fig_rho_orig)

    # Interpolated spatial random effects
    fig_rho = far.plot.rho("output_roadless/rho.tif",
                           borders="data/ctry_PROJ.shp",
                           output_file="output_roadless/rho.png")
    plt.close(fig_rho)

    # Spatial probability of deforestation
    fig_prob = far.plot.prob("output_roadless/prob.tif",
                             borders="data/ctry_PROJ.shp",
                             output_file="output_roadless/prob.png")
    plt.close(fig_prob)

    # Forest-cover change 2015-2050
    fig_fcc_35yr = far.plot.fcc("output_roadless/fcc_35yr.tif",
                                borders="data/ctry_PROJ.shp",
                                output_file="output_roadless/fcc_35yr.png")
    plt.close(fig_fcc_35yr)
    
    # Forest-cover change 2015-2025
    fig_fcc_10yr = far.plot.fcc("output_roadless/fcc_10yr.tif",
                                borders="data/ctry_PROJ.shp",
                                output_file="output_roadless/fcc_10yr.png")
    plt.close(fig_fcc_10yr)
    
    # Forest-cover change 2015-2050 with half deforestation
    fig_fcc_35yr_half = far.plot.fcc("output_roadless/fcc_35yr_half.tif",
                                     borders="data/ctry_PROJ.shp",
                                     output_file="output_roadless/fcc_35yr_half.png")
    plt.close(fig_fcc_35yr_half)

# End
