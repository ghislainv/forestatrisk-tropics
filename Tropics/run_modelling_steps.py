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
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import log_loss
import pandas as pd
import re  # Regular expressions


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

    :param fcc_source: Source of forest cover data. Either "jrc" or
    "gfc".

    """


    # Make output directory
    far.make_dir("output")

    # ========================================================
    # Sample points
    # ========================================================

    dataset = far.sample(nsamp=10000, adapt=True,
                         Seed=1234, csize=10,
                         var_dir="data",
                         input_forest_raster="fcc23.tif",
                         output_file="output/sample.txt",
                         blk_rows=0)

    # # Import data as pandas DataFrame if necessary
    # dataset = pd.read_table("output/sample.txt", delimiter=",")

    # Remove NA from data-set (otherwise scale() and
    # model_binomial_iCAR doesn't work)
    dataset = dataset.dropna(axis=0)
    # Set number of trials to one for far.model_binomial_iCAR()
    dataset["trial"] = 1

    # Sample size
    ndefor = sum(dataset.fcc23 == 0)
    nfor = sum(dataset.fcc23 == 1)
    with open("output/sample_size.csv", "w") as f:
        f.write("var, n\n")
        f.write("ndefor, " + str(ndefor) + "\n")
        f.write("nfor, " + str(nfor) + "\n")

    # Descriptive statistics
    # Model formulas
    formula_1 = "fcc23 ~ dist_road + dist_town + dist_river + \
    dist_defor + dist_edge + altitude + slope - 1"
    # Standardized variables (mean=0, std=1)
    formula_2 = "fcc23 ~ scale(dist_road) + scale(dist_town) + \
    scale(dist_river) + scale(dist_defor) + scale(dist_edge) + \
    scale(altitude) + scale(slope) - 1"
    formulas = (formula_1, formula_2)

    # Exceptions
    if iso3 == "VIR":  # No river for VIR
        # Model formula
        formula_1 = "fcc23 ~ dist_road + dist_town + \
        dist_defor + dist_edge + altitude + slope - 1"
        # Standardized variables (mean=0, std=1)
        formula_2 = "fcc23 ~ scale(dist_road) + scale(dist_town) + \
        scale(dist_defor) + scale(dist_edge) + \
        scale(altitude) + scale(slope) - 1"
        formulas = (formula_1, formula_2)

    # Loop on formulas
    for f in range(len(formulas)):
        # Output file
        of = "output/correlation_" + str(f) + ".pdf"
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
    with open("output/summary_hSDM.txt", "w") as f:
        f.write(str(mod_binomial_iCAR))

    # Plot
    figs = mod_binomial_iCAR.plot(output_file="output/mcmc.pdf",
                                  plots_per_page=3,
                                  figsize=(9, 6),
                                  dpi=300)
    plt.close("all")

    # ========================================================
    # Model performance comparison: cross-validation
    # ========================================================

    # Cross-validation for icar, glm and RF
    CV_df_icar = far.cross_validation(dataset, formula, mod_type="icar", ratio=30, nrep=5,
                                      icar_args={"n_neighbors": nneigh, "neighbors": adj,
                                                 "burnin": 1000, "mcmc": 1000, "thin": 1,
                                                 "beta_start": mod_binomial_iCAR.betas})
    CV_df_glm = far.cross_validation(dataset, formula, mod_type="glm", ratio=30, nrep=5)
    CV_df_rf = far.cross_validation(dataset, formula, mod_type="rf", ratio=30, nrep=5,
                                    rf_args={"n_estimators": 500, "n_jobs": 3})

    # Save result to disk
    CV_df_icar.to_csv("output/CV_icar.csv", header=True, index=False)
    CV_df_glm.to_csv("output/CV_glm.csv", header=True, index=False)
    CV_df_rf.to_csv("output/CV_rf.csv", header=True, index=False)
    
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

    # Random forest model
    formula_rf = formula
    y, x = dmatrices(formula_rf, data=dataset, NA_action="drop")
    Y = y[:, 0]
    X_rf = x[:, :-1]  # We remove the last column (cells)
    mod_rf = RandomForestClassifier(n_estimators=500,
                                    n_jobs=3)
    mod_rf = mod_rf.fit(X_rf, Y)

    # Deviances
    deviance_null = 2*log_loss(Y, mod_null.predict_proba(X_null), normalize=False)
    deviance_glm = 2*log_loss(Y, mod_glm.predict_proba(X_glm), normalize=False)
    deviance_rf = 2*log_loss(Y, mod_rf.predict_proba(X_rf), normalize=False)
    deviance_icar = mod_binomial_iCAR.deviance
    deviance_full = 0
    dev = [deviance_null, deviance_glm, deviance_rf, deviance_icar, deviance_full]

    # Result table
    mod_dev = pd.DataFrame({"model": ["null", "glm", "rf", "icar", "full"],
                            "deviance": dev})
    perc = 100*(1-mod_dev.deviance/deviance_null)
    mod_dev["perc"] = perc
    mod_dev = mod_dev.round(0)
    mod_dev.to_csv("output/model_deviance.csv", header=True, index=False)

    # ========================================================
    # Interpolating spatial random effects
    # ========================================================

    # Spatial random effects
    rho = mod_binomial_iCAR.rho

    # Interpolate
    far.interpolate_rho(rho=rho, input_raster="data/fcc23.tif",
                        output_file="output/rho.tif",
                        csize_orig=10, csize_new=1)

    # ========================================================
    # Predicting spatial probability of deforestation
    # ========================================================

    # Update dist_edge and dist_defor at t3
    os.rename("data/dist_edge.tif", "data/dist_edge.tif.bak")
    os.rename("data/dist_defor.tif", "data/dist_defor.tif.bak")
    copy2("data/forecast/dist_edge_forecast.tif", "data/dist_edge.tif")
    copy2("data/forecast/dist_defor_forecast.tif", "data/dist_defor.tif")

    # Compute predictions
    far.predict_raster_binomial_iCAR(
        mod_binomial_iCAR, var_dir="data",
        input_cell_raster="output/rho.tif",
        input_forest_raster="data/forest/forest_t3.tif",
        output_file="output/prob.tif",
        blk_rows=10  # Reduced number of lines to avoid memory problems
    )

    # Reinitialize data
    os.remove("data/dist_edge.tif")
    os.remove("data/dist_defor.tif")
    os.rename("data/dist_edge.tif.bak", "data/dist_edge.tif")
    os.rename("data/dist_defor.tif.bak", "data/dist_defor.tif")

    # ========================================================
    # Mean annual deforestation rate (ha.yr-1)
    # ========================================================

    # Forest cover
    p = re.compile("BRA-.*")
    m = p.match(iso3)
    if m is not None:  # For Brazil
        if (fcc_source == "jrc"):
            fcc_BRA = pd.read_csv("../fcc_BRA_jrc.csv")
        else:
            fcc_BRA = pd.read_csv("../fcc_BRA_gfc.csv")
    else:  # Other country
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
        # Annual deforestation
        T = 10.0
        annual_defor = (fc[2] - fc[4]) / T

    # Dates and time intervals
    dates_fut = ["2030", "2035", "2040", "2050", "2055", "2060", "2070", "2080", "2085", "2090", "2100"]
    ndates_fut = len(dates_fut)
    ti = [10, 15, 20, 30, 35, 40, 50, 60, 65, 70, 80]

    # ========================================================
    # Predicting forest cover change
    # ========================================================

    # Loop on dates
    for i in range(ndates_fut):
        # Amount of deforestation (ha)
        if m is not None:  # For Brazil
            defor = fcc_BRA.loc[fcc_BRA["iso3"]==iso3, "defor" + dates_fut[i]].values[0]
        else:
            defor = np.rint(annual_defor * ti[i])
        # Compute future forest cover
        stats = far.deforest(input_raster="output/prob.tif",
                             hectares=defor,
                             output_file="output/fcc_" + dates_fut[i] + ".tif",
                             blk_rows=128)
        # Save some stats if date = 2050
        if dates_fut[i] == "2050":
            # Save stats to disk with pickle
            pickle.dump(stats, open("output/stats.pickle", "wb"))
            # Plot histograms of probabilities
            fig_freq = far.plot.freq_prob(stats,
                                          output_file="output/freq_prob.png")
            plt.close(fig_freq)

    # ========================================================
    # Carbon emissions
    # ========================================================

    # Create dataframe
    dpast = ["2020"]
    dpast.extend(dates_fut)
    C_df = pd.DataFrame({"date": dpast, "C": np.repeat(-99, ndates_fut + 1)},
                        columns=["date","C"])
    # Loop on date
    for i in range(ndates_fut):
        carbon = far.emissions(input_stocks="data/emissions/AGB.tif",
                               input_forest="output/fcc_" + dates_fut[i] + ".tif")
        C_df.loc[C_df["date"]==dates_fut[i], ["C"]] = carbon
    # Past emissions 
    carbon = far.emissions(input_stocks="data/emissions/AGB.tif",
                           input_forest="data/fcc23.tif")
    C_df.loc[C_df["date"]==dpast[0], ["C"]] = carbon
    # Save dataframe
    C_df.to_csv("output/C_emissions.csv", header=True, index=False)

    # ========================================================
    # Figures
    # ========================================================

    # Forest at t3
    fig_forest = far.plot.forest("data/forest/forest_t3.tif",
                                 maxpixels=1e8,
                                 borders="data/ctry_PROJ.shp",
                                 output_file="output/forest_t3.png")
    plt.close(fig_forest)

    # Forest-cover change 123
    fig_fcc = far.plot.fcc123("data/forest/fcc123.tif",
                              maxpixels=1e8,
                              borders="data/ctry_PROJ.shp",
                              output_file="output/fcc123.png")
    plt.close(fig_fcc)

    # Forest-cover change 12345
    fig_fcc = far.plot.fcc12345("data/forest/fcc12345.tif",
                                maxpixels=1e8,
                                borders="data/ctry_PROJ.shp",
                                output_file="output/fcc12345.png")
    plt.close(fig_fcc)

    # Original spatial random effects
    fig_rho_orig = far.plot.rho("output/rho_orig.tif",
                                borders="data/ctry_PROJ.shp",
                                output_file="output/rho_orig.png")
    plt.close(fig_rho_orig)

    # Interpolated spatial random effects
    fig_rho = far.plot.rho("output/rho.tif",
                           borders="data/ctry_PROJ.shp",
                           output_file="output/rho.png")
    plt.close(fig_rho)

    # Spatial probability of deforestation
    fig_prob = far.plot.prob("output/prob.tif",
                             maxpixels=1e8,
                             borders="data/ctry_PROJ.shp",
                             output_file="output/prob.png")
    plt.close(fig_prob)

    # Projected future forest-cover change
    for i in range(ndates_fut):
        fig_fcc = far.plot.fcc("output/fcc_" + dates_fut[i] + ".tif",
                               maxpixels=1e8,
                               borders="data/ctry_PROJ.shp",
                               output_file="output/fcc_" + dates_fut[i] + ".png")
        plt.close(fig_fcc)

# End
