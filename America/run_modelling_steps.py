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
def run_modelling_steps(iso3, fcc_source="jrc"):
    """Runs the modelling steps.

    Runs the modelling steps: (1) select the significant variables,
    (2) model the deforestation process, (3) interpolate the spatial
    random effects, (3) predict the probability of deforestation, (4)
    compute mean annul deforested areas, and (4) forecast the future
    forest cover at different date in the future following a
    business-as-usual scenario.

    :param iso3: Country ISO 3166-1 alpha-3 code. This is used to
    handle exceptions in the modelling process for specific countries.

    :param fcc_source: Source of the forest cover data to compute
    time-interval for the observation of the past
    deforestation. Either "jrc" or "gfc".

    """


    # Make output directory
    far.make_dir("output_jrc")

    # ========================================================
    # Sample points
    # ========================================================
    
    # dataset = far.sample(nsamp=10000, adapt=True,
    #                      Seed=1234, csize=10,
    #                      var_dir="data",
    #                      input_forest_raster="fcc23.tif",
    #                      output_file="output_jrc/sample.txt",
    #                      blk_rows=0)

    # Import data as pandas DataFrame if necessary
    dataset = pd.read_table("output_jrc/sample.txt", delimiter=",")

    # Remove NA from data-set (otherwise scale() and
    # model_binomial_iCAR doesn't work)
    dataset = dataset.dropna(axis=0)
    # Set number of trials to one for far.model_binomial_iCAR()
    dataset["trial"] = 1

    # Sample size
    ndefor = sum(dataset.fcc23 == 0)
    nfor = sum(dataset.fcc23 == 1)
    with open("output_jrc/sample_size.csv", "w") as f:
        f.write("var, n\n")
        f.write("ndefor, " + str(ndefor) + "\n")
        f.write("nfor, " + str(nfor) + "\n")

    # # Descriptive statistics
    # # Model formulas
    # formula_1 = "fcc23 ~ dist_road + dist_town + dist_river + \
    # dist_defor + dist_edge + altitude + slope - 1"
    # # Standardized variables (mean=0, std=1)
    # formula_2 = "fcc23 ~ scale(dist_road) + scale(dist_town) + \
    # scale(dist_river) + scale(dist_defor) + scale(dist_edge) + \
    # scale(altitude) + scale(slope) - 1"
    # formulas = (formula_1, formula_2)

    # # Exceptions
    # if iso3 == "VIR":  # No river for VIR
    #     # Model formula
    #     formula_1 = "fcc23 ~ dist_road + dist_town + \
    #     dist_defor + dist_edge + altitude + slope - 1"
    #     # Standardized variables (mean=0, std=1)
    #     formula_2 = "fcc23 ~ scale(dist_road) + scale(dist_town) + \
    #     scale(dist_defor) + scale(dist_edge) + \
    #     scale(altitude) + scale(slope) - 1"
    #     formulas = (formula_1, formula_2)

    # # Loop on formulas
    # for f in range(len(formulas)):
    #     # Output file
    #     of = "output_jrc/correlation_" + str(f) + ".pdf"
    #     # Data
    #     y, data = dmatrices(formulas[f], data=dataset,
    #                         return_type="dataframe")
    #     # Plots
    #     figs = far.plot.correlation(y=y, data=data,
    #                                 plots_per_page=3,
    #                                 figsize=(7, 8),
    #                                 dpi=300,
    #                                 output_file=of)
    #     plt.close("all")

    # ========================================================
    # hSDM model
    # ========================================================

    # Spatial cells for spatial-autocorrelation
    nneigh, adj = far.cellneigh(raster="data/fcc23.tif", csize=10, rank=1)

    # List of variables
    variables = ["C(pa)", "scale(altitude)", "scale(slope)",
                 "scale(dist_defor)", "scale(dist_edge)", "scale(dist_road)",
                 "scale(dist_town)", "scale(dist_river)"]
    # Exceptions
    if iso3 == "VIR": variables.remove("scale(dist_river)")
    if iso3 == "SXM": variables.remove("C(pa)")
    # Transform into numpy array
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
    with open("output_jrc/summary_hSDM.txt", "w") as f:
        f.write(str(mod_binomial_iCAR))

    # Plot
    figs = mod_binomial_iCAR.plot(output_file="output_jrc/mcmc.pdf",
                                  plots_per_page=3,
                                  figsize=(9, 6),
                                  dpi=300)
    plt.close("all")

    # ========================================================
    # Model performance comparison: cross-validation
    # ========================================================

    # Cross-validation for icar and glm
    CV_df_icar = far.cross_validation(dataset, formula, mod_type="icar", ratio=30, nrep=5,
                                      icar_args={"n_neighbors": nneigh, "neighbors": adj,
                                                 "burnin": 1000, "mcmc": 1000, "thin": 1,
                                                 "beta_start": mod_binomial_iCAR.betas})
    CV_df_glm = far.cross_validation(dataset, formula, mod_type="glm", ratio=30, nrep=5)

    # Save result to disk
    CV_df_icar.to_csv("output_jrc/CV_icar.csv", header=True, index=False)
    CV_df_glm.to_csv("output_jrc/CV_glm.csv", header=True, index=False)
    
    # ========================================================
    # Model performance comparison: deviance
    # ========================================================

    # Null model
    formula_null = "I(1-fcc23) ~ 1"
    y, x = dmatrices(formula_null, data=dataset, NA_action="drop")
    Y = y[:, 0]
    X_null = x[:, :]
    mod_null = LogisticRegression(solver="lbfgs")
    mod_null = mod_null.fit(X_null, Y)

    # Simple glm with no spatial random effects
    formula_glm = formula
    y, x = dmatrices(formula_glm, data=dataset, NA_action="drop")
    Y = y[:, 0]
    X_glm = x[:, :-1]  # We remove the last column (cells)
    mod_glm = LogisticRegression(solver="lbfgs")
    mod_glm = mod_glm.fit(X_glm, Y)

    # Deviances
    deviance_null = 2*log_loss(Y, mod_null.predict_proba(X_null), normalize=False)
    deviance_glm = 2*log_loss(Y, mod_glm.predict_proba(X_glm), normalize=False)
    deviance_icar = mod_binomial_iCAR.deviance
    deviance_full = 0
    dev = [deviance_null, deviance_glm, deviance_icar, deviance_full]

    # Result table
    mod_dev = pd.DataFrame({"model": ["null", "glm", "icar", "full"],
                            "deviance": dev})
    perc = 100*(1-mod_dev.deviance/deviance_null)
    mod_dev["perc"] = perc
    mod_dev = mod_dev.round(0)
    mod_dev.to_csv("output_jrc/model_deviance.csv", header=True, index=False)

    # # ========================================================
    # # Interpolating spatial random effects
    # # ========================================================

    # # Spatial random effects
    # rho = mod_binomial_iCAR.rho

    # # Interpolate
    # far.interpolate_rho(rho=rho, input_raster="data/fcc23.tif",
    #                     output_file="output_jrc/rho.tif",
    #                     csize_orig=10, csize_new=1)

    # # ========================================================
    # # Predicting spatial probability of deforestation
    # # ========================================================

    # # Update dist_edge and dist_defor at t3
    # os.rename("data/dist_edge.tif", "data/dist_edge.tif.bak")
    # os.rename("data/dist_defor.tif", "data/dist_defor.tif.bak")
    # copy2("data/forecast/dist_edge_forecast.tif", "data/dist_edge.tif")
    # copy2("data/forecast/dist_defor_forecast.tif", "data/dist_defor.tif")

    # # Compute predictions
    # far.predict_raster_binomial_iCAR(
    #     mod_binomial_iCAR, var_dir="data",
    #     input_cell_raster="output_jrc/rho.tif",
    #     input_forest_raster="data/forest/forest_t3.tif",
    #     output_file="output_jrc/prob.tif",
    #     blk_rows=128
    # )

    # # Reinitialize data
    # os.remove("data/dist_edge.tif")
    # os.remove("data/dist_defor.tif")
    # os.rename("data/dist_edge.tif.bak", "data/dist_edge.tif")
    # os.rename("data/dist_defor.tif.bak", "data/dist_defor.tif")

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
                                 maxpixels=1e8,
                                 borders="data/ctry_PROJ.shp",
                                 output_file="output_jrc/forest_t3.png")
    plt.close(fig_forest)

    # Forest-cover change 2000-2019
    fig_fcc = far.plot.fcc("data/forest/fcc13.tif",
                           maxpixels=1e8,
                           borders="data/ctry_PROJ.shp",
                           output_file="output_jrc/fcc13.png")
    plt.close(fig_fcc)

    # Forest-cover change 2010-2019
    fig_fcc = far.plot.fcc("data/fcc23.tif",
                           maxpixels=1e8,
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
                             maxpixels=1e8,
                             borders="data/ctry_PROJ.shp",
                             output_file="output_jrc/prob.png")
    plt.close(fig_prob)

    # Projected future forest-cover change
    for i in range(ndate):
        fig_fcc = far.plot.fcc("output_jrc/fcc_" + date[i] + ".tif",
                               maxpixels=1e8,
                               borders="data/ctry_PROJ.shp",
                               output_file="output_jrc/fcc_" + date[i] + ".png")
        plt.close(fig_fcc)

# End
