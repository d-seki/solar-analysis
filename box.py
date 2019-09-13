import matplotlib.pyplot as plt
import numpy as np
from sunpy.coordinates.utils import GreatArc
from astropy.coordinates import SkyCoord
import astropy.units as u


def map_box(amap, bl_arc, tr_arc, img=None, c='c'):
    color = c
    if type(bl_arc) != SkyCoord:
        bl_arc = SkyCoord(
            bl_arc[0]*u.arc, bl_arc[1]*u.arc,
            frame=amap.coordinate_frame)
        tr_arc = SkyCoord(
            tr_arc[0]*u.arc, tr_arc[1]*u.arc,
            frame=amap.coordinate_frame)
    br_arc = SkyCoord(
	tr_arc.Tx, bl_arc.Ty,
	frame=amap.coordinate_frame)
    tl_arc = SkyCoord(
	bl_arc.Tx, tr_arc.Ty,
	frame=amap.coordinate_frame)
    blbr =  GreatArc(bl_arc, br_arc)
    brtr =  GreatArc(br_arc, tr_arc)
    tltr =  GreatArc(tl_arc, tr_arc)
    bltl =  GreatArc(bl_arc, tl_arc)
    if not img:
        img = amap.plot()

    img.axes.plot_coord(blbr.coordinates(), color=color)
    img.axes.plot_coord(brtr.coordinates(), color=color)
    img.axes.plot_coord(tltr.coordinates(), color=color)
    img.axes.plot_coord(bltl.coordinates(), color=color)

    return img


def plot_box(img, bl, tr, c='black'):
    # bl[x,y], tr[x,y]
    img.axes.hlines(bl[1], bl[0],tr[0], color=c)
    img.axes.hlines(tr[1], bl[0],tr[0], color=c)
    img.axes.vlines(bl[0], bl[1],tr[1], color=c)
    img.axes.vlines(tr[0], bl[1],tr[1], color=c)

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

