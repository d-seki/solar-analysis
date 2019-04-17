# coding: utf-8
# version: Python 3.6
import numpy as np

def MapArcToPix(hmap, xra, yra):
    return np.array(
            [xra / hmap.scale[0].value + hmap.reference_pixel[0].value,\
            yra / hmap.scale[1].value + hmap.reference_pixel[1].value])

def FitsArcToPix(tmpfits, xra, yra):
    return np.array(
            [xra/tmpfits[1].header['cdelt1']+tmpfits[1].header['crpix1'] ,\
            yra/tmpfits[1].header['cdelt2']+tmpfits[1].header['crpix2']])

