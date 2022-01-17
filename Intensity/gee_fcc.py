#!/usr/bin/python
# -*- coding: utf-8 -*-

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ghislainv.github.io
# python_version  :>=2.7
# license         :GPLv3
# ==============================================================================

# PRIOR TO EXECUTING THE FOLLOWING SCRIPT
# AUTHENTICATE TO GOOGLE EARTH ENGINE
# Use following command line: earthengine authenticate

# Annual product legend
# ---------------------
# 1. Undisturbed Tropical moist forest (TMF)
# 2. Degraded TMF
# 3. Deforested land
# 4. Forest regrowth
# 5. Permanent or seasonal water
# 6. Other land cover

# Imports
import ee
# import geemap

# Initialize
ee.Initialize()

continent = ["America", "Africa", "Asia"]
cont = ["Ame", "Afr", "Asi"]
ncont = len(continent)

# ===============================
# Function to compute forest area
# ===============================


# Compute_forest_area
def compute_forest_area(feature):
    vals = forest.multiply(ee.Image.pixelArea()).reduceRegion(
        reducer=ee.Reducer.sum(), geometry=feature.geometry(),
        scale=30, maxPixels=5e9)
    return feature.set(vals).select([".*"], None, False)


# Loop on continents
for c in range(ncont):

    # =========
    # Data
    # =========

    # JRC annual product (AP)
    AP = ee.ImageCollection("projects/JRC/TMF/v1_2020/AnnualChanges")
    AP = AP.mosaic().toByte()

    # Country grid
    ctry_grid = ee.FeatureCollection("users/ghislainv/"
                                     "forestatrisk-tropics-intensity/"
                                     "borders_" + continent[c] + "_grid")

    # ==============================
    # Forest cover from 2000 to 2020
    # ==============================

    # ap_allYear: forest if Y = 1 or 2.
    AP_forest = AP.where(AP.eq(2), 1)
    ap_allYear = AP_forest.where(AP_forest.neq(1), 0)

    # Forest cover Jan 2000
    ap_2000_2021 = ap_allYear.select(list(range(9, 31)))
    forest2000 = ap_2000_2021.reduce(ee.Reducer.sum()).gte(1)
    forest = forest2000.rename('fc2000')

    # Loop on year (max Janv 2021)
    for i in range(1, 22):
        ap_20XX_2021 = ap_allYear.select(list(range(9 + i, 31)))
        forest20XX = ap_20XX_2021.reduce(ee.Reducer.sum()).gte(1)
        band = forest20XX.rename('fc' + str(2000 + i))
        forest = ee.Image.cat(forest, band)

    # ===============================
    # Extract some countries
    # ===============================

    # # extractCountry
    # def extractCountry(area_code):
    #     return ctry_grid.filterMetadata("area_code", "equals", area_code)

    # countryList = ee.List(["MDG", "GNQ"])
    # some_countries = ee.FeatureCollection(
    #     countryList.map(extractCountry)
    # ).flatten()

    # Country grid for computations (calculations)
    # ctry_grid_calc = some_countries
    # ctry_grid_calc = ctry_grid

    # Map the counry borders with grid
    # Map = geemap.Map(center=[40, -100], zoom=4)
    # Map.addLayer(ctry_grid_calc)
    # Map

    # ===================================
    # Maps an algorithm over a collection
    # ===================================

    # Compute forest area
    forest_area = ctry_grid.map(compute_forest_area)

    # ===================================
    # Export
    # ===================================

    # Export table to drive
    task = ee.batch.Export.table.toDrive(
        collection=forest_area,
        description='ForestArea_grid_' + cont[c],
        folder='forestatrisk-tropics-intensity',
        fileNamePrefix='ForestArea_grid_' + cont[c])
    task.start()

# EOF
