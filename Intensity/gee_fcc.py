#!/usr/bin/python
# -*- coding: utf-8 -*-

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ghislainv.github.io
# python_version  :2.7
# license         :GPLv3
# ==============================================================================

# PRIOR TO EXECUTING THE FOLLOWING SCRIPT
# AUTHENTICATE TO GOOGLE EARTH ENGINE
# Use following command line: earthengine authenticate

# Annual product legend
# ---------------------
#  1. Tropical moist forest (TMF including bamboo-dominated forest and
#     mangroves)
#  2. TMF converted later in a tree plantation
#  3. NEW degradation
#  4. Ongoing degradation (disturbances still detected - can be few years after
#     first degrad if several degradation stages)
#  5. Degraded forest (former degradation, no disturbances detected anymore)
#  6. NEW deforestation (may follow degradation)
#  7. Ongoing deforestation (disturbances still detected)
#  8. NEW Regrowth
#  9. Regrowthing
# 10. Other land cover (not water)
# 11. Permanent Water (pekel et al.2015)
# 12. Seasonal Water (pekel et al.2015)
# 13. Init period without valid data - Init class = TMF
# 14. Init period with min 1 valid obs - Init class = TMF
# 15. Nodata  - Init class = other LC
# 16. Init period without valid data - Init class = Plantation

# Imports
import ee
# import geemap

# Initialize
ee.Initialize()

# =========
# Data
# =========

# JRC annual product (AP)
# AP = ee.ImageCollection(path + "AnnualChanges1982_2019")
AP = ee.ImageCollection("users/ClassifLandsat072015/"
                        "Roadless2019/AnnualChanges_1982_2019")
AP = AP.mosaic().toByte()

# Country grid
ctry_grid = ee.FeatureCollection("users/ghislainv/"
                                 "forestatrisk-tropics-intensity/"
                                 "borders_Asia_grid")

# ==============================
# Forest cover from 2000 to 2020
# ==============================

# ap_allYear: forest if Y = 1, 2, 3, 4, 5, 13, or 14.
AP_forest = AP.where(AP.eq(2).Or(AP.eq(3)).Or(AP.eq(4)).Or(
    AP.eq(5)).Or(AP.eq(13)).Or(AP.eq(14)), 1)
ap_allYear = AP_forest.where(AP_forest.neq(1), 0)

# Forest cover Jan 2000
ap_2000_2020 = ap_allYear.select(list(range(17, 38)))
forest2000 = ap_2000_2020.reduce(ee.Reducer.sum()).gte(1)
forest = forest2000.rename('fc2000')

# Loop on year
for i in range(1, 21):
    ap_20XX_2020 = ap_allYear.select(list(range(17 + i, 38)))
    forest20XX = ap_20XX_2020.reduce(ee.Reducer.sum()).gte(1)
    band = forest20XX.rename('fc' + str(2000 + i))
    forest = ee.Image.cat(forest, band)

# ===============================
# Function to compute forest area
# ===============================


# compute_forest_area
def compute_forest_area(feature):
    vals = forest.multiply(ee.Image.pixelArea()).reduceRegion(
        reducer=ee.Reducer.sum(), geometry=feature.geometry(),
        crs="EPSG:3395", scale=30, maxPixels=5e9)
    return feature.set(vals).select([".*"], None, False)


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
# Maps an algorithm over a collection
# ===================================
# Export table to drive
task = ee.batch.Export.table.toDrive(
    collection=forest_area,
    description='ForestArea_grid_Asi',
    folder='forestatrisk-tropics-intensity',
    fileNamePrefix='ForestArea_grid_Asi')
task.start()

# EOF
