forestatrisk-tropics
================

[![Website
ForestAtRisk](https://img.shields.io/badge/website-ForestAtRisk-blue)](https://forestatrisk.cirad.fr)
[![License
GPLv3](https://img.shields.io/badge/licence-GPLv3-8f10cb.svg)](https://www.gnu.org/licenses/gpl-3.0.html)

This repository includes the code used to produce the results of the
following scientific article:

<a href="https://orcid.org/0000-0002-1685-4997"><img alt="ORCID logo" src="Website/images/Logo_ORCID.png" width="16" height="16" /></a>
<a href="https://ecology.ghislainv.fr" style="color:#2C3E50;">**Vieilledent
G.,**</a>
<a href="https://orcid.org/0000-0003-3851-8588"><img alt="ORCID logo" src="Website/images/Logo_ORCID.png" width="16" height="16" /></a>
<a href="https://www.researchgate.net/profile/Christelle_Vancutsem" style="color:#2C3E50;">**C.
Vancutsem,**</a>
<a href="https://www.researchgate.net/profile/Achard_Frederic" style="color:#2C3E50;">**F.
Achard.**</a> Spatial forecasting of forest cover change in the humid
tropics over the 21<sup>st</sup> century. in prep.

## Reproducibility of the results

### Computing historical deforestation rates and uncertainty

``` bash
Rscript Intensity/intensity.R
```

### Geoprocessing and modelling

1\. Derive historical forest cover change map from the annual product of
Vancutsem et al. 2021 using Google Earth Engine.

``` bash
python combine.py
```

2\. Download raw data from on-line databases (GADM, SRTM, WDPA, OSM),
and Google Drive.

``` bash
combine.py
```

3\. Compute explanatory variables (elevation, slope, distances, etc.).

``` bash
combine.py
```

4\. Model and forecast

``` bash
combine.py
```

5\. Combine rasters to obtain continental maps

``` bash
Maps/combine.py
```

### R for post-processing and writing

1\. Plot main maps

``` bash
Rscript Maps/main_maps.R
```

2\. Plot supplementary maps

``` bash
Rscript Maps/supp_maps.R
```

3\. Synthesize results

``` bash
Rscript Analysis/synthesis.R
```

4\. Compile documents

``` bash
Rscript Manuscript/zzz_knitr_compile/compile_book.R
```

<!-- End of file -->
