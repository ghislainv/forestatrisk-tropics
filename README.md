
[![License
GPLv3](https://img.shields.io/badge/licence-GPLv3-8f10cb.svg)](https://www.gnu.org/licenses/gpl-3.0.html)
[![Cirad
Dataverse](https://img.shields.io/badge/DOI-10.18167/DVN1/7N2BTU-green)](https://doi.org/10.18167/DVN1/7N2BTU)
[![Website
ForestAtRisk](https://img.shields.io/badge/web-ForestAtRisk-blue)](https://forestatrisk.cirad.fr)
[![forestatrisk Python
package](https://img.shields.io/badge/python-forestatrisk-306998?logo=python&logoColor=ffd43b&color=306998)](https://ecology.ghislainv.fr/forestatrisk)

This repository includes the code used to produce the results of the
following scientific article:

<a href="https://orcid.org/0000-0002-1685-4997"><img alt="ORCID logo" src="Website/images/Logo_ORCID.png" width="16" height="16" /></a>
**Vieilledent G.,**
<a href="https://orcid.org/0000-0003-3851-8588"><img alt="ORCID logo" src="Website/images/Logo_ORCID.png" width="16" height="16" /></a>
**C. Vancutsem,** **F. Achard.** Spatial forecasting of forest cover
change in the humid tropics over the 21<sup>st</sup> century. in prep.

<img alt="ORCID logo" src="Manuscript/Article/figures/prob_zoom.png" />

Figure 2: **Pantropical map of the risk of deforestation.**

## Minimal reproducible example using the `forestatrisk` Python package

This
[notebook](https://ecology.ghislainv.fr/forestatrisk/notebooks/far_tropics.html)
provides a minimal and reproducible example presenting the general
approach we followed to model and forecast deforestation in each of the
119 study areas (representing 92 countries) considered in the above
article. We use the Guadeloupe archipelago as a case study. The notebook
is available at the [website](https://ecology.ghislainv.fr/forestatrisk)
associated with the `forestatrisk` Python package. This package has been
specifically developed for this study and provides functions to model
and forecast deforestation in the tropics.

## Steps followed to produce the results of the study

We present below the R and Python scripts which have been used to
produce the results of the study, from the datasets preparation to the
writing of the manuscript.

### 1. Preparing datasets

``` bash
## Derive past forest cover change maps from the annual product 
## of Vancutsem et al. 2021 using Google Earth Engine.
python Tropics/forest_gee_jrc.py

## Download raw data from on-line databases (GADM, SRTM, WDPA, OSM), and Google Drive.
python Tropics/download_raw_data.py

## Compute explanatory variables (elevation, slope, distances, etc.).
python Tropics/compute_variables.py
```

### 2. Estimating deforestation intensity

``` bash
## Compute deforestation rates and uncertainty
Rscript Intensity/intensity.R

## Estimate contagious deforestation between states of Brazil
python Intensity/brazil_fcc_jrc.py
```

### 3. Spatial modeling and forecasting

``` bash
## Model and forecast
python Tropics/model_and_forecast.py
```

### 4. Post-processing and writing

``` bash
## Combine rasters to obtain continental maps
python Maps/combine.py

## Plot main maps
Rscript Maps/main_maps.R

## Plot supplementary maps
Rscript Maps/supp_maps.R

## Synthesize results
Rscript Analysis/synthesis.R

## Compile documents
Rscript Manuscript/zzz_knitr_compile/compile_book.R
```

## Website accompanying the article

A website at <https://forestatrisk.cirad.fr> is accompanying the above
article. The website includes the following ressources:

### Interactive map

Interactive maps from this study (forest cover change, deforestation
risk, and projected forest cover in 2050 and 2100) have been made
available:

-   [Map of the tropics](https://forestatrisk.cirad.fr/maps.html)

### Download

Rasters of results from this study can be downloaded as Cloud Optimized
GeoTIFFs ([COG](https://www.cogeo.org/)):

-   [Rasters](https://forestatrisk.cirad.fr/rasters.html)
-   [COG tutorial](https://forestatrisk.cirad.fr/notebooks/cog.html)

### Supplementary data

-   [Data S1](https://forestatrisk.cirad.fr/data-s.html): Uncertainty
    around projected forest cover by study-area.
-   [Data S2](https://forestatrisk.cirad.fr/data-s.html): Uncertainty
    around cumulative carbon emissions associated with future
    deforestation.

### `forestatrisk` Python package

Results from this study have been obtained with the `forestatrisk`
Python package:

-   [Package website](https://ecology.ghislainv.fr/forestatrisk/) (with
    full documentation)
-   [Tutorials](https://ecology.ghislainv.fr/forestatrisk/articles.html)

<span style="display: block; height: 15px;"></span>
<p>
Copyright © 2021 <a href="https://www.cirad.fr/en/">Cirad</a>,
<a href="https://ec.europa.eu/jrc/en">EC JRC</a>. All rights reserved.
</p>

<a href="https://www.cirad.fr/en/"><img alt="RF" src="Website/images/Logo_RF.jpg" height="75"></a>
<a href="https://www.cirad.fr/en/"><img alt="Cirad" src="Website/images/Logo_Cirad.jpg" height="60"></a>
<a href="https://amap.cirad.fr"><img alt="AMAP" src="Website/images/Logo_AMAP.jpg" height="60"></a>
<a href="https://ec.europa.eu/jrc/en"><img alt="Cirad" src="Website/images/Logo_JRC.png" height="60"></a>

<!-- End of file -->
