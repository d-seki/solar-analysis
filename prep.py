# -*- coding: utf-8 -*-
# version: Python 3.6

import astropy.units as u
import numpy as np
from sunpy.map import Map
from sunpy.io import fits
from disk import disk
from sunpy import coordinates as coord
import drms


def url_aiaprep(utime,wavelength,imgsize=2048, ratio=1.2, margin=100):
    wavelength = str(wavelength)
    c = drms.Client(email='seki@kusastro.kyoto-u.ac.jp', verbose=True)
    aiaquery = c.query(
        'aia.lev1_euv_12s['+utime.strftime('%Y.%m.%d_%H:%M:%S')+']'\
        '[?wavelnth='+wavelength+'?]', seg='image')
    url = 'http://jsoc.stanford.edu' + aiaquery.image[0]
    newmap = aiaprep(Map(url),
                imgsize=imgsize,ratio=ratio,margin=margin)
    return newmap

def url_hmiprep(utime, imgsize=2048, ratio=1.2, radius=0.95, padding=True):
    c = drms.Client(email='seki@kusastro.kyoto-u.ac.jp', verbose=True)
    hmiquery = c.query(
        'hmi.M_45s['+utime.strftime('%Y.%m.%d_%H:%M:%S_TAI')+']',
        seg='magnetogram', key='**ALL**')
    url = 'http://jsoc.stanford.edu' + hmiquery[1]
    dum = Map(url.values[0,0])
    dum.meta.update(hmiquery[0].T.to_dict()[0])
    newmap = hmiprep(dum,
                imgsize=imgsize,ratio=ratio,radius=radius,padding=padding)
    newmap.meta.update({
        'date-obs': newmap.meta['date__obs']
        })
    return newmap

def synptcprep(synptcmap):
    synptcmap.meta['CUNIT2'] = 'degree'
    synptcmap.meta['CUNIT1'] = 'degree'
    synptcmap.meta['CDELT2'] = 180./np.pi*synptcmap.meta['CDELT2']
    synptcmap.meta['CDELT1'] *= -1
    synptcmap.meta.update({'date-obs':synptcmap.meta['date']})
    return synptcmap

def aiaprep(aiamap, imgsize=2048, ratio=1.2, margin=100):
    newdata = np.zeros((imgsize, imgsize), 'f')
    centpixsize = int(imgsize/2)
    expcrrctdata = aiamap.data/aiamap.meta['exptime']
    expcrrctmap = Map(expcrrctdata, aiamap.meta)
    rawsolradpix = aiamap.meta['rsun_obs']/aiamap.scale[0].value
    solradpix = imgsize/2/ratio
    rotmap = expcrrctmap.rotate(angle=aiamap.meta['crota2']*u.degree,
            scale=solradpix/rawsolradpix, recenter=True, missing=0)
    centerpix = [round(rotmap.reference_pixel[0].value),
                round(rotmap.reference_pixel[1].value)]
    newdata[
            centpixsize-int(round(solradpix+margin)):\
            centpixsize+int(round(solradpix+margin)),
            centpixsize-int(round(solradpix+margin)):\
            centpixsize+int(round(solradpix+margin))
            ] = rotmap.data[
                    centerpix[0]-int(round(solradpix+margin)):\
                    centerpix[0]+int(round(solradpix+margin)),
                    centerpix[1]-int(round(solradpix+margin)):\
                    centerpix[1]+int(round(solradpix+margin))
                    ]
    newmap = Map(newdata, rotmap.meta)

    newmap.meta['crpix1'] = int(imgsize/2)
    newmap.meta['crpix2'] = int(imgsize/2)
    newmap.meta['exptime'] = 1.
    newmap.meta['lvl_num'] = 1.5

    return newmap


def hmiprep(hmimap, imgsize=2048, ratio=1.2, radius=0.95, margin=0, padding=True):
    newdata = np.zeros((imgsize, imgsize), 'f')
    centpixsize = int(imgsize/2)
    rawsolradpix = hmimap.meta['rsun_obs']/hmimap.scale[0].value
    solradpix = imgsize/2/ratio
    mask = disk([imgsize, imgsize], solradpix*radius)
    rotmap = hmimap.rotate(angle=hmimap.meta['crota2']*u.degree,
            scale=solradpix/rawsolradpix, recenter=True, missing=0)
    centerpix = [round(rotmap.reference_pixel[0].value),
                round(rotmap.reference_pixel[1].value)]
    newdata[
            centpixsize-int(round(solradpix+margin)):\
            centpixsize+int(round(solradpix+margin)),
            centpixsize-int(round(solradpix+margin)):\
            centpixsize+int(round(solradpix+margin))
            ] = rotmap.data[
                    centerpix[0]-int(round(solradpix+margin)):\
                    centerpix[0]+int(round(solradpix+margin)),
                    centerpix[1]-int(round(solradpix+margin)):\
                    centerpix[1]+int(round(solradpix+margin))
                    ]
    newmap = Map(newdata*mask, rotmap.meta)
    if not padding:
        mask0 = disk([imgsize, imgsize], solradpix)
        newmap.data[np.where(mask0 == 0)] = np.nan

    newmap.meta['crpix1'] = int(imgsize/2)
    newmap.meta['crpix2'] = int(imgsize/2)
    newmap.meta['lvl_num'] = 1.5

    return newmap


def smartprep(hafile, imgsize=2048, ratio=1.2, margin=50):
    hadata = fits.read(hafile)[0]
    hadata.header['DATE'] = hadata.header['DATE'].split('Z')[0]
    hadata.header['DATE-OBS'] = hadata.header['DATE-OBS'].split('Z')[0]
    hadata.header.update({
        'crota2': hadata.header['CROTA'],
        'radius': hadata.header['R0'],
        'wavelnth': hadata.header['WAVE_LEN'],
        'dsun_obs': coord.get_sunearth_distance(
                    hadata.header['DATE-OBS']).to(u.m).value,
        'crln_obs': hadata.header['L0'],
        'crlt_obs': hadata.header['B0']*180/np.pi,
        'cunit1':'arcsec',
        'cunit2':'arcsec',
        'ctype1': 'HPLN-TAN',
        'ctype2': 'HPLT-TAN',
        'X0': hadata.header['XCEN'],
        'Y0': hadata.header['YCEN'],
        })

    hamap = Map([hadata.data, hadata.header])
    newdata = np.zeros((imgsize, imgsize), 'f')
    centpixsize = int(imgsize/2)
    rawsolradpix = hamap.meta['r0']/hamap.scale[0].value
    solradpix = imgsize/2/ratio
    rotmap = hamap.rotate(angle=hamap.meta['crota']*u.degree,
            scale=solradpix/rawsolradpix, recenter=True, missing=0)
    centerpix = [round(rotmap.reference_pixel[0].value),
                round(rotmap.reference_pixel[1].value)]
    newdata[
            centpixsize-int(round(solradpix+margin)):\
            centpixsize+int(round(solradpix+margin)),
            centpixsize-int(round(solradpix+margin)):\
            centpixsize+int(round(solradpix+margin))
            ] = rotmap.data[
                    centerpix[0]-int(round(solradpix+margin)):\
                    centerpix[0]+int(round(solradpix+margin)),
                    centerpix[1]-int(round(solradpix+margin)):\
                    centerpix[1]+int(round(solradpix+margin))
                    ]
    newmap = Map(newdata, rotmap.meta)
    newmap.meta['crpix1'] = int(imgsize/2)
    newmap.meta['crpix2'] = int(imgsize/2)
    newmap.meta['lvl_num'] = 1.5

    return newmap

# deprecated
#def smartprep(hamap, imgsize=2048, ratio=1.2, margin=50):
#    hamap.meta['date'] = hamap.meta['date'].split('Z')[0]
#    hamap.meta['date-obs'] = hamap.meta['date-obs'].split('Z')[0]
#    hamap.meta.update({
#        'crota2': hamap.meta['crota'],
#        'radius': hamap.meta['r0'],
#        'wavelnth': hamap.meta['wave_len'],
#        'dsun_obs': coord.get_sunearth_distance(
#                    hamap.meta['date-obs']).to(u.m).value,
#        'crln_obs': hamap.meta['l0'],
#        'crlt_obs': hamap.meta['b0']*180/np.pi,
#        'cunit1':'arcsec',
#        'cunit2':'arcsec',
#        'ctype1': 'HPLN-TAN',
#        'ctype2': 'HPLT-TAN',
#        'X0': hamap.meta['xcen'],
#        'Y0': hamap.meta['ycen'],
#        })
#
#    newdata = np.zeros((imgsize, imgsize), 'f')
#    centpixsize = int(imgsize/2)
#    rawsolradpix = hamap.meta['r0']/hamap.scale[0].value
#    solradpix = imgsize/2/ratio
#    rotmap = hamap.rotate(angle=hamap.meta['crota']*u.degree,
#            scale=solradpix/rawsolradpix, recenter=True, missing=0)
#    centerpix = [round(rotmap.reference_pixel[0].value),
#                round(rotmap.reference_pixel[1].value)]
#    newdata[
#            centpixsize-int(round(solradpix+margin)):\
#            centpixsize+int(round(solradpix+margin)),
#            centpixsize-int(round(solradpix+margin)):\
#            centpixsize+int(round(solradpix+margin))
#            ] = rotmap.data[
#                    centerpix[0]-int(round(solradpix+margin)):\
#                    centerpix[0]+int(round(solradpix+margin)),
#                    centerpix[1]-int(round(solradpix+margin)):\
#                    centerpix[1]+int(round(solradpix+margin))
#                    ]
#    newmap = Map(newdata, rotmap.meta)
#    newmap.meta['crpix1'] = int(imgsize/2)
#    newmap.meta['crpix2'] = int(imgsize/2)
#    newmap.meta['lvl_num'] = 1.5
#
#    return newmap


def smartprep2(hafits, num, imgsize=2048, ratio=1.2):
    """
    hafits = fits.open(file)
    """
    hamap = Map(hafits[0].data[num], hafits[0].header)
    hamap.meta['date-obs'] = hamap.meta['date-obs'].split('Z')[0]
    hamap.meta.update({
        'date': hamap.meta['date-obs'],
        'crota2': hamap.meta['crota'],
        'radius': hamap.meta['r0'],
        'wavelnth': hamap.meta['wave_len'],
        'dsun_obs': coord.get_sunearth_distance(
                    hamap.meta['date-obs']).to(u.m).value,
        'crln_obs': hamap.meta['l0'],
        'crlt_obs': hamap.meta['b0']*180/np.pi
        })
    newdata = np.zeros((imgsize, imgsize), 'f')
    centpixsize = int(imgsize/2)
    rawsolradpix = hamap.meta['r0']/hamap.scale[0].value
    solradpix = imgsize/2/ratio
    rotmap = hamap.rotate(angle=hamap.meta['crota']*u.degree,
            scale=solradpix/rawsolradpix, recenter=True, missing=0)
    centerpix = [round(rotmap.reference_pixel[0].value), round(rotmap.reference_pixel[1].value)]
    newdata[
            centpixsize-round(solradpix):centpixsize+round(solradpix),
            centpixsize-round(solradpix):centpixsize+round(solradpix)
            ] = rotmap.data[
            centerpix[0]-round(solradpix):centerpix[0]+round(solradpix),
            centerpix[1]-round(solradpix):centerpix[1]+round(solradpix)
            ]
    newmap = Map(newdata, rotmap.meta)
    newmap.meta['crpix1'] = int(imgsize/2)
    newmap.meta['crpix2'] = int(imgsize/2)
    newmap.meta['lvl_num'] = 1.5

    return newmap


def gongprep(hamap, imgsize=2048, ratio=1.2):
    hamap.meta.update({
        'crota2': hamap.meta['roll'],
        'radius': hamap.meta['solar-r'],
        'dsun_obs': coord.get_sunearth_distance(
                    hamap.meta['date-obs']).to(u.m).value,
        'crln_obs': hamap.meta['solar-l0'],
        'crlt_obs': hamap.meta['solar-b0'], # degree ??
        'cdelt1': hamap.meta['solar-r']/hamap.meta['radius'],
        'cdelt2': hamap.meta['solar-r']/hamap.meta['radius']
        })
    hamap.meta['ctype1'] = 'HPLN-TAN'
    hamap.meta['ctype2'] = 'HPLT-TAN'

    newdata = np.zeros((imgsize, imgsize), 'f')
    centpixsize = int(imgsize/2)
    rawsolradpix = hamap.meta['radius']/hamap.scale[0].value
    solradpix = imgsize/2/ratio
    rotmap = hamap.rotate(angle=hamap.meta['crota2']*u.degree,
            scale=solradpix/rawsolradpix, recenter=True, missing=0)
    centerpix = [round(rotmap.reference_pixel[0].value),
                round(rotmap.reference_pixel[1].value)]
    newdata[
            centpixsize-round(solradpix):centpixsize+round(solradpix),
            centpixsize-round(solradpix):centpixsize+round(solradpix)
            ] = rotmap.data[
            centerpix[0]-round(solradpix):centerpix[0]+round(solradpix),
            centerpix[1]-round(solradpix):centerpix[1]+round(solradpix)
            ]
    newmap = Map(newdata, rotmap.meta)
    newmap.meta['crpix1'] = int(imgsize/2)
    newmap.meta['crpix2'] = int(imgsize/2)
    newmap.meta['lvl_num'] = 1.5

    return newmap

def norhprep(norhmap, imgsize=512, ratio=1.2):
    norhmap.meta.update({
        'wavelnth': norhmap.meta['obs-freq'],
        'date-obs': norhmap.meta['date-obs']+' '+norhmap.meta['time-obs'],
        'crota2': 0,
        'dsun_obs': coord.get_sunearth_distance(
                    norhmap.meta['date-obs']).to(u.m).value,
        'crlt_obs': norhmap.meta['solb'],
        'radius': norhmap.meta['solr']
        })
    norhmap.meta['ctype1'] = 'HPLN-TAN'
    norhmap.meta['ctype2'] = 'HPLT-TAN'

    newdata = np.zeros((imgsize, imgsize), 'f')
    centpixsize = int(imgsize/2)
    rawsolradpix = norhmap.meta['radius']/norhmap.scale[0].value
    solradpix = imgsize/2/ratio
    rotmap = norhmap.rotate(angle=norhmap.meta['crota2']*u.degree,
            scale=solradpix/rawsolradpix, recenter=True, missing=0)
    centerpix = [round(rotmap.reference_pixel[0].value),
            round(rotmap.reference_pixel[1].value)]
    newdata[
            centpixsize-round(solradpix):centpixsize+round(solradpix),
            centpixsize-round(solradpix):centpixsize+round(solradpix)
            ] = rotmap.data[
            centerpix[0]-round(solradpix):centerpix[0]+round(solradpix),
            centerpix[1]-round(solradpix):centerpix[1]+round(solradpix)
            ]
    newmap = Map(newdata, rotmap.meta)
    newmap.meta['crpix1'] = int(imgsize/2)
    newmap.meta['crpix2'] = int(imgsize/2)
    newmap.meta['lvl_num'] = 1.5

    return newmap






"""

**** sunpy.map.Mapは、まだまだ発展途上！えらくつかいづらいので注意！ ****

cdelt1,2 --> x-y軸のarcsec表示に寄与。pixelに対し、何arcsecか。太陽部分の画像は変わらない。
太陽部分の画像を変えたければ、rotateでscaleを変えるしかない。

"""
