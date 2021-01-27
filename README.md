forestatrisk-tropics
================

This repository is accompanying the following scientific article:

**Vieilledent G., C. Vancutsem, F. Achard.** Spatial forecasting of
forest cover change in the humid tropics over the 21<sup>st</sup>
century. in prep.

## Reproducibility of the results

### Python for geoprocessing and modelling

1\. Derive historical forest cover change map from the annual product of
Vancutsem et al. 2020 on Google Earth Engine.

``` bash
combine.py
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

6\. Plot main maps

``` bash
Rscript Maps/main_maps.R
```

7\. Plot supplementary maps

``` bash
Rscript Maps/supp_maps.R
```

8\. Synthesize results

``` bash
Rscript Analysis/synthesis.R
```

9\. Compile documents

``` bash
Rscript Manuscript/zzz_knitr_compile/compile_book.R
```

## Associated ressources

### Interactive maps

Interactive maps from this study (forest cover change, deforestation
risk, and projected forest cover in 2050 and 2100) have been made
available:

-   [Tropics](https://forestatrisk.cirad.fr/tropics/)
-   [Madagascar](https://forestatrisk.cirad.fr/mada/)
-   [New-Caledonia](https://forestatrisk.cirad.fr/newcal/)

### Downloading

Rasters of results from this study can be downloaded as Cloud Optimized
GeoTIFF ([COG](https://www.cogeo.org/)):

-   [Rasters](https://forestatrisk.cirad.fr/rasters.html)
-   [COG tutorial](https://forestatrisk.cirad.fr/notebooks/cog.html)

### Supplementary data

-   [Data S1](https://forestatrisk.cirad.fr/tropics/supplementary-data):
    Past and projected forest cover.
-   [Data S2](https://forestatrisk.cirad.fr/tropics/supplementary-data):
    Cumulative carbon emissions associated to future deforestation

### `forestatrisk` Python package

Results from this study have been obtained with the `forestatrisk`
Python package:

-   [Package website](https://ecology.ghislainv.fr/forestatrisk/) (with
    full documentation)
-   [Tutorials](https://ecology.ghislainv.fr/forestatrisk/articles.html)

<!-- End of file -->
