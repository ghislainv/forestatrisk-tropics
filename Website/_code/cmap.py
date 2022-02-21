#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ecology.ghislainv.fr
# python_version  :>=3
# license         :GPLv3
# ==============================================================================

# This script create a color map for terracotta
# see: https://terracotta-python.readthedocs.io/en/latest/tutorials/custom-colormaps.html

# Standard library imports
import os
import subprocess

# Third party imports
import numpy as np


# Rescale function
def rescaling(x, range1, range2):
    """Rescale a value from range1 to range2"""
    delta1 = range1[1] - range1[0]
    delta2 = range2[1] - range2[0]
    val = (delta2 * (x - range1[0]) / delta1) + range2[0]
    return np.int_(val)


# Color interpolation
def colorFader(c1, c2, mix=0):
    """Fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)"""
    c1 = np.array(c1)
    c2 = np.array(c2)
    return (1-mix)*c1 + mix*c2


# Set wd
os.chdir("/home/ghislain/Code/forestatrisk-tropics/Website/_code")

# Deforestation breaks in the interval 1-255
rval = rescaling(np.array([1, 39322, 52429, 65535]),
                 [1, 65535], [1, 255])
# >>> rval
# array([  1, 153, 204, 255])

# Number of colors between breaks (including breaks)
ncol = (rval[1:] - rval[:-1]) + 1
# >>> ncol
# array([153,  52,  52])

# RGB colors for breaks
green = (34, 139, 34, 255)
orange = (255, 165, 0, 255)
red = (227, 26, 28, 255)
black = (0, 0, 0, 255)

# Array for colors (must be type unint8 for terracotta)
cmap_data = np.zeros(shape=(255, 4), dtype="uint8")

# Loops on the number of colors
for n in range(ncol[0]):
    cmap_data[n] = colorFader(green, orange, mix=n/(ncol[0] - 1))

for n in range(ncol[1]):
    offset = ncol[0] - 1
    cmap_data[n + offset] = colorFader(orange, red, mix=n/(ncol[1] - 1))

for n in range(ncol[2]):
    offset = ncol[0] + ncol[1] - 2
    cmap_data[n + offset] = colorFader(red, black, mix=n/(ncol[2] - 1))

# Save color maps
np.save('prob_rgba.npy', cmap_data)

# Copy to fdb
fdb_path = ("fdb:/home/www/forestatrisk/tropics/"
            "tc-cmaps/prob_rgba.npy")
cmd = f"scp prob_rgba.npy {fdb_path}"
subprocess.call(cmd, shell=True)


# Convert to hex color for Google Earth Engine
def hextriplet(colortuple):
    return '#' + ''.join(f'{i:02X}' for i in colortuple)


hex_list = [hextriplet(x) for x in cmap_data]
print(hex_list)

# EOF
