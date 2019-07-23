#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from glob import glob
import matplotlib.pyplot as plt
import os
import matplotlib.colors as colors
import cartopy.crs as ccrs


datadir = os.path.join(os.path.split(__file__)[0], 'data')
figdir = os.path.join(os.path.split(__file__)[0], 'figs')

#files = glob(os.path.join(datadir, '*s0012*fits'))
#files = glob(os.path.join(datadir, '*4-4*fits'))
files = glob(os.path.join(datadir, '*fits'))
files.sort()

# central longitude of the projection
cenlon = 0.
# if the projection has a dividing line where we have to split
mustsplit = False #89.5 W, 66.2 S
tr = ccrs.Orthographic(central_longitude=-89.5, central_latitude=-66.2)
#tr = ccrs.Mollweide()

# minimum and maximum flux for the colorbar
vmin = 100
vmax = 1001.
cnorm = colors.LogNorm(vmin=vmin, vmax=vmax)
cmap = 'gray'

test = True
makefig = True
doclean = False
highres = False
savefig = False

fname = 'ortho.png'
savefile = os.path.join(figdir, fname)


credit = 'By Ethan Kruse\n@ethan_kruse'
##################################################################

# XXX: 
#files[0] = files[-2]
#files = files[-3:]


def grab_sector(sector):
    file = os.path.join(datadir, 'tesscurl_sector_{0}_ffic.sh'.format(sector))
    with open(file, 'r') as ff:
        lines = ff.readlines()
    lines.sort()
    
    midorbit = len(lines)*7//8
    ll = lines[midorbit]
    date = ll.split('tess')[1].split('-')[0]
    ret = [line for line in lines if date in line]
    assert len(ret) == 16
    for iret in ret:
        print(iret)
    return


def clean(data):
    import scipy.ndimage
    fdata = data * 1
    fdata[fdata > 1000] = 100
    fdata[fdata < 0] = 0
    fdata = scipy.ndimage.uniform_filter(fdata, size=201)
    # fdata[fdata > 300] = 300
    
    """
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(lat[1:,1:], lon[1:,1:], fdata, cmap='gray',
                       linewidth=0, antialiased=False)
    """
    return data - fdata + 100


plt.close('all')

for ii, ifile in enumerate(files):
    with fits.open(ifile) as ff:
        wcs = WCS(ff[1].header)
            
        data = ff[1].data * 1

        xinds = np.arange(-0.5, data.shape[0]-0.4)
        yinds = np.arange(-0.5, data.shape[1]-0.4)
        mesh = np.meshgrid(xinds, yinds, indexing='ij')
        
        lon, lat = wcs.all_pix2world(mesh[1].flatten(), mesh[0].flatten(), 0)
        
        lon = lon.reshape(mesh[0].shape)
        lat = lat.reshape(mesh[1].shape)

        # chop out unexposed rows/columns
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
        
        if doclean:
            data = clean(data)
        
        print(f'{ii+1} of {len(files)}: {lon.min():.2f}, {lon.max():.2f}, {lat.min():.2f}, {lat.max():.2f}')
        
        if makefig:
            if ii == 0 or test:
                if highres:
                    fig = plt.figure(figsize=(80,80))
                else:
                    fig = plt.figure()
                ax = plt.axes([0.01, 0.01, 0.99, 0.99], projection=tr)
                data_tr = ccrs.PlateCarree()
                # ax.coastlines()
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
            if highres:
                fsz = 150
            else:
                fsz = 14
            plt.text(0.02, 0.02, credit, transform=fig.transFigure, ha='left', 
                     va='bottom', multialignment='left', fontsize=fsz, fontname='Carlito')

if makefig and savefig:
    inum = 1
    orig = savefile
    while os.path.exists(savefile):
        savefile = os.path.splitext(orig)[0] + f'{inum}' + os.path.splitext(orig)[1]
        inum += 1
    plt.savefig(savefile)

"""
from mpl_toolkits.mplot3d import Axes3D
import scipy.ndimage

fdata = data * 1
fdata[fdata > 1000] = 100
fdata = scipy.ndimage.uniform_filter(fdata, size=201)


fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(lat[1:,1:], lon[1:,1:], fdata, cmap='gray',
                       linewidth=0, antialiased=False)



plt.figure()
ax = plt.axes([0.01, 0.01, 0.99, 0.99], projection=tr)
plt.pcolormesh(lon, lat, data, norm=cnorm, alpha=0.3, transform=data_tr, cmap=cmap)
"""



