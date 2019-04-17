import matplotlib.pyplot as plt
import numpy as np
from sunpy.coordinates.utils import GreatArc
from astropy.coordinates import SkyCoord
import astropy.units as u


def map_box(amap,bl_arcsec, tr_arcsec, ax=None, color='c'):
    bl_arcsec = SkyCoord(
        bl_arcsec[0]*u.arcsec, bl_arcsec[1]*u.arcsec,
	frame=amap.coordinate_frame)
    tr_arcsec = SkyCoord(
	tr_arcsec[0]*u.arcsec, tr_arcsec[1]*u.arcsec,
	frame=amap.coordinate_frame)
    br_arcsec = SkyCoord(
	tr_arcsec.Tx, bl_arcsec.Ty,
	frame=amap.coordinate_frame)
    tl_arcsec = SkyCoord(
	bl_arcsec.Tx, tr_arcsec.Ty,
	frame=amap.coordinate_frame)
    blbr =  GreatArc(bl_arcsec, br_arcsec)
    brtr =  GreatArc(br_arcsec, tr_arcsec)
    tltr =  GreatArc(tl_arcsec, tr_arcsec)
    bltl =  GreatArc(bl_arcsec, tl_arcsec)
    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    amap.plot()
    ax.plot_coord(blbr.coordinates(), color=color)
    ax.plot_coord(brtr.coordinates(), color=color)
    ax.plot_coord(tltr.coordinates(), color=color)
    ax.plot_coord(bltl.coordinates(), color=color)


def plot_box(img, bl, tr, clr='black'):
    bl = bl[::-1]
    tr = tr[::-1]
    plt.hlines(bl[0], bl[1],tr[1], color=clr)
    plt.hlines(tr[0], bl[1],tr[1], color=clr)
    plt.vlines(bl[1], bl[0],tr[0], color=clr)
    plt.vlines(tr[1], bl[0],tr[0], color=clr)

def draw_box(img, bl, tr, clr='white', num=10):
    width = int(img.shape[0]/200)
    if clr == 'black':
        result = img.copy()
        result[bl[0]:tr[0],bl[1]-2:bl[1]+2] = -1
        result[bl[0]:tr[0],tr[1]-2:tr[1]+2] = -1
        result[bl[0]-2:bl[0]+2,bl[1]:tr[1]] = -1
        result[tr[0]-2:tr[0]+2,bl[1]:tr[1]] = -1

    if clr == 'white':
        result = img.copy()
        Ly = int((tr[0] - bl[0])/num)
        Lx = int((tr[1] - bl[1])/num)
        for i in range(0,num,2):
            result[bl[0]+Ly*i:bl[0]+Ly*(i+1),bl[1]-width:bl[1]+width] = np.nan
            result[bl[0]+Ly*i:bl[0]+Ly*(i+1),tr[1]-width:tr[1]+width] = np.nan
        for j in range(0,num,2):
            result[bl[0]-width:bl[0]+width,bl[1]+Lx*j:bl[1]+Lx*(j+1)] = np.nan
            result[tr[0]-width:tr[0]+width,bl[1]+Lx*j:bl[1]+Lx*(j+1)] = np.nan

    return result

