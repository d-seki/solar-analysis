# -*- coding: utf-8 -*-

import numpy as np
import math

def disk(window, radius):
    y,x = np.ogrid[0:window[0],0:window[1]]
    disk = np.sqrt(((y-window[0]/2)**2 + (x-window[1]/2)**2)) <= radius

    return disk
