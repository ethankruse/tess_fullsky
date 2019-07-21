#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 11:20:39 2019

@author: ekruse1
"""

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from glob import glob
import matplotlib.pyplot as plt
import os
import matplotlib.colors as colors
from astropy.modeling import projections
import cartopy.crs as ccrs

datadir = '/Users/ekruse1/research/tess_sky/data'

# files = glob(os.path.join(datadir, '*s0001*fits'))
files = glob(os.path.join(datadir, '*fits'))
files.sort()
# minimum and maximum flux for the regular movies
vmin = 70
vmax = 1001.

# central longitude of the projection
cenlon = 0.
# if the projection has a dividing line where we have to split
mustsplit = False

cnorm = colors.LogNorm(vmin=vmin, vmax=vmax)
cmap = 'gray'


def grab_sector(sector):
    file = '/Users/ekruse1/research/tess_sky/data/tesscurl_sector_{0}_ffic.sh'.format(sector)
    with open(file, 'r') as ff:
        lines = ff.readlines()
    lines.sort()
    
    midorbit = len(lines)//4
    ll = lines[midorbit]
    date = ll.split('tess')[1].split('-')[0]
    ret = [line for line in lines if date in line]
    assert len(ret) == 16
    return ret


#lines = grab_sector(1)
#for line in lines:
#    print(line)


plt.close('all')

#plt.figure()

for ii, ifile in enumerate(files):
    with fits.open(ifile) as ff:
        wcs = WCS(ff[1].header)
        
        if ii == 0:
            plt.figure()
            tr = ccrs.Orthographic(central_longitude=179., central_latitude=-70.)
            #tr = ccrs.Mollweide()
            ax = plt.axes(projection=tr)
            data_tr = ccrs.PlateCarree()
            ax.coastlines()
            
        data = ff[1].data * 1
        # print(ff[0].header['date-obs'])

        xinds = np.arange(-0.5, data.shape[0]-0.4)
        yinds = np.arange(-0.5, data.shape[1]-0.4)
        mesh = np.meshgrid(xinds, yinds, indexing='ij')
        
        lon, lat = wcs.all_pix2world(mesh[1].flatten(), mesh[0].flatten(), 0)
        
        lon = lon.reshape(mesh[0].shape)
        lat = lat.reshape(mesh[1].shape)

        goodx = np.where(np.median(data, axis=1) > 50.)[0]
        # make sure we're not cutting out rows in the middle
        goodx = np.arange(goodx[0], goodx[-1]+0.5).astype(int)
        data = data[goodx, :]
        uind = np.unique(np.concatenate((goodx, goodx+1)))
        lon = lon[uind, :]
        lat = lat[uind, :]
        
        goody = np.where(np.median(data, axis=0) > 50.)[0]
        goody = np.arange(goody[0], goody[-1]+0.5).astype(int)
        data = data[:, goody]
        uind2 = np.unique(np.concatenate((goody, goody+1)))
        lon = lon[:, uind2]
        lat = lat[:, uind2]
        # lon must be between -180 and 180 instead of 0 to 360
        lon -= 180.
        
        print(ii, lon.min(), lon.max(), lat.min(), lat.max())
        # for wraparounds:
        if lon.max() > cenlon + 120 and lon.min() < cenlon - 120 and mustsplit:
            left = np.where(lon > cenlon + 120)
            lonleft = lon * 1
            lonleft[left] = cenlon - 180.
            plt.pcolormesh(lonleft, lat, data, norm=cnorm, alpha=0.3, transform=data_tr, cmap=cmap)
            
            right = np.where(lon < cenlon - 120)
            lonright = lon * 1
            lonright[right] = cenlon + 180.
            plt.pcolormesh(lonright, lat, data, norm=cnorm, alpha=0.3, transform=data_tr, cmap=cmap)
        else:
            plt.pcolormesh(lon, lat, data, norm=cnorm, alpha=0.3, transform=data_tr, cmap=cmap)
        #plt.text(np.median(lon), np.median(lat), '{0}'.format(ii), transform=data_tr)

        