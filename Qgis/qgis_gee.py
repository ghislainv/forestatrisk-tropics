#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ecology.ghislainv.fr
# python_version  :>=2.7
# license         :GPLv3
# ==============================================================================

# Annual product legend
# 1. Tropical moist forest (TMF including bamboo-dominated forest and mangroves)
# 2. TMF converted later in a tree plantation
# 3. NEW degradation 
# 4. Ongoing degradation (disturbances still detected - can be few years after first degrad if several degradation stages)
# 5. Degraded forest (former degradation, no disturbances detected anymore)
# 6. NEW deforestation (may follow degradation)
# 7. Ongoing deforestation (disturbances still detected)
# 8. NEW Regrowth
# 9. Regrowthing
# 10. Other land cover (not water)
# 11. Permanent Water (pekel et al.2015)     
# 12. Seasonal Water (pekel et al.2015) 
# 13. Init period without valid data - Init class = TMF
# 14. Init period with min 1 valid obs - Init class = TMF
# 15. Nodata  - Init class = other LC 
# 16. Init period without valid data - Init class = Plantation

# Imports
import ee
from ee_plugin import Map

# rgb2hex
def rgb2hex(r, g, b):
    hex = "#{:02x}{:02x}{:02x}".format(r, g, b)
    return hex

# Extent
iso3 = "REU"
#iso3 = "MDG"
if iso3 == "REU":
    extent_latlong = [55.216251, -21.389860, 55.837360, -20.871805]
elif iso3 == "MDG":
    extent_latlong = [43.10, -25.15, 51.70, -11.80]
elif iso3 == "BEN":
    extent_latlong = [0.7745, 6.2351, 3.8517, 12.4183]
else:
    extent_latlong = [49.307035, -16.028757, 50.526517, -14.810834]  # Masoala

# Region
region = ee.Geometry.Rectangle(extent_latlong, proj="EPSG:4326",
                               geodesic=False)
region = region.buffer(10000).bounds()
export_coord = region.getInfo()["coordinates"]

# Path to JRC products
#path = "users/ghislainv/jrc/"

# JRC annual product (AP)
#AP = ee.ImageCollection(path + "AnnualChanges1982_2019")
AP = ee.ImageCollection(
    "users/ClassifLandsat072015/Roadless2019/AnnualChanges_1982_2019")
AP = AP.mosaic().toByte().unmask()
AP = AP.clip(region)
MASK = AP.select(0)
AP = AP.where(MASK.eq(0), 0)

## temp
PALETTEAnnualChange = [
rgb2hex(0,0,0), 
rgb2hex(0,90,0),     # val 1. Tropical moist forest (TMF including bamboo-dominated forest and mangroves)
rgb2hex(0,110,0),    # val 2. TMF converted later in a tree plantation
rgb2hex(155,80,60),     # val 3. NEW degradation 
rgb2hex(135,115,45),        # val 4. Ongoing degradation (disturbances still detected - can be few years after first degrad if several degradation stages)
rgb2hex(100,135,35),        # val 5. Degraded forest (former degradation, no disturbances detected anymore)
rgb2hex(255,20,0),      # val 6. NEW deforestation (may follow degradation)
rgb2hex(255,255,155),       # val 7. Ongoing deforestation (disturbances still detected)
rgb2hex(152,230,0),     # val 8. NEW Regrowth
rgb2hex(50,160,0),          # val 9. Regrowthing
rgb2hex(255,255,255), # val 10. Other land cover (not water)
rgb2hex(0,77,168),    # val 11. Permanent Water (pekel et al.2015)     
rgb2hex(0,157,200),   # val 12. Seasonal Water (pekel et al.2015) 
rgb2hex(0,80,0),      # val 13. Init period without valid data - Init class = TMF
rgb2hex(0,80,0),      # val 14. Init period with min 1 valid obs - Init class = TMF
rgb2hex(255,255,255), # val 15. Nodata  - Init class = other LC 
rgb2hex(0,80,0),      # val 16. Init period without valid data - Init class = Plantation
]

## Map.addLayer(AP.select('Dec2010'),{'min':0, 'max': 16, 'palette': PALETTEAnnualChange}, "Statuf of one year AnnualChanges 2010", True)

## /temp

# ap_allYear: forest if Y = 1, 2, 3, 4, 5, 13, or 14.
AP_forest = AP.where(AP.eq(2).Or(AP.eq(3)).Or(AP.eq(4)).Or(
    AP.eq(5)).Or(AP.eq(13)).Or(AP.eq(14)), 1)
ap_allYear = AP_forest.where(AP_forest.neq(1), 0)

# Forest in Jan 2020
forest2020 = ap_allYear.select(37)

# Map.addLayer(forest2020, {'min': 0, 'max': 1, 'palette': ["red", "green"]}, "forest2020", True)

# Forest cover Jan 2015
ap_2015_2020 = ap_allYear.select(list(range(32, 38)))
forest2015 = ap_2015_2020.reduce(ee.Reducer.sum())
forest2015 = forest2015.gte(1)

# Forest cover Jan 2010
ap_2010_2020 = ap_allYear.select(list(range(27, 38)))
forest2010 = ap_2010_2020.reduce(ee.Reducer.sum())
forest2010 = forest2010.gte(1)

# Forest cover Jan 2005
ap_2005_2020 = ap_allYear.select(list(range(22, 38)))
forest2005 = ap_2005_2020.reduce(ee.Reducer.sum())
forest2005 = forest2005.gte(1)

# Forest cover Jan 2000
ap_2000_2020 = ap_allYear.select(list(range(17, 38)))
forest2000 = ap_2000_2020.reduce(ee.Reducer.sum())
forest2000 = forest2000.gte(1)

# Forest raster with five bands
forest = forest2000.addBands(forest2005).addBands(
    forest2010).addBands(forest2015).addBands(forest2020)
forest = forest.select([0, 1, 2, 3, 4], ["forest2000", "forest2005",
                                         "forest2010", "forest2015",
                                         "forest2020"])
forest = forest.set("system:bandNames", ["forest2000", "forest2005",
                                         "forest2010", "forest2015",
                                         "forest2020"])

# Forest raster with three values
forest = forest2000.where(forest2010.eq(1), 2).where(forest2020.eq(1), 3)

# Palette
Palette = [rgb2hex(255, 255, 255), rgb2hex(255, 165, 0), rgb2hex(227, 26, 28), rgb2hex(34, 139, 34)]

# Map
Map.addLayer(forest.updateMask(forest), {'min': 0, 'max': 3, "palette": Palette}, "forest2000-2020", True)
if iso3 == "REU":
    x, y, z = (55.53, -21.11, 12)
elif iso3 == "MDG":
    x, y, z = (46.70, -19.30, 8)
elif iso3 == "BEN":
    x, y, z = (0.7745, 6.2351, 8)
else:
    x, y, z = (50.30, -15.47, 12)
Map.setCenter(x, y, z)

# EOF