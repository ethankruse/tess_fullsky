#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from astropy.coordinates import SkyCoord, BarycentricTrueEcliptic
import pandas as pd
from astropy import units as u

cornerdir = os.path.join(os.path.split(__file__)[0], 'corners')
figdir = os.path.join(os.path.split(__file__)[0], 'figs')

useexoplots = True
if useexoplots:
    toilist = os.path.join(os.path.split(__file__)[0], 'exoplots_data.csv')
else:
    toilist = os.path.join(os.path.split(__file__)[0], 'toi_list.txt')
# options are 'north', 'south', or 'both'
hemisphere = 'north'
# for full-sky Mollweide projections, do we want to use ecliptic coordinates
# if False, uses celestial coordinates (ICRS, right ascension/declination)
ecliptic_coords = True
# make a Plate Caree image instead of Mollweide
platecarree = False

# for hemisphere == 'both', whether to include the Kepler/K2 footprints
addkepler = True

# option to not print any text on the images
notext = False

# parameters that change depending on the hemisphere
if hemisphere == 'both':
    # coordinates at the center of the projection
    cenlon = 0.
    cenlat = 0.
    # if the projection has a dividing line where we have to wrap data around
    wrap = True
    # set up our desired map projection
    if platecarree:
        tr = ccrs.PlateCarree(central_longitude=cenlon)
        fbase = 'platecarree'
    else:
        tr = ccrs.Mollweide(central_longitude=cenlon)
        # what the output file name base should be
        fbase = 'mollweide'
    # the coordinates of the corners of the CCDs
    edgefiles = [os.path.join(os.path.split(__file__)[0], 'edges_south.txt'),
                 os.path.join(os.path.split(__file__)[0], 'edges_north.txt')]
    if ecliptic_coords:
        fbase += '_ecliptic'
    if addkepler:
        fbase += '_kepler'
    #  title text in upper left
    if notext:
        title = ''
    else:
        title = "NASA TESS's View\nof the Sky"
elif hemisphere == 'south':
    # south ecliptic pole coordinates are 90, -66.560708333333
    cenlon = 90.
    cenlat = -66.560708333333
    wrap = False
    tr = ccrs.AzimuthalEquidistant(central_longitude=cenlon,
                                   central_latitude=cenlat)
    edgefiles = [os.path.join(os.path.split(__file__)[0], 'edges_south.txt')]
    fbase = 'azeq_south'
    if notext:
        title = ''
    else:
        title = "NASA TESS's View\nof the Southern\nHemisphere"
    # turn off ecliptic coordinates since it doesn't matter
    ecliptic_coords = False
    # for now, don't try to put Kepler/K2 on the hemisphere maps
    addkepler = False
elif hemisphere == 'north':
    cenlon = -90.
    cenlat = 66.560708333333
    wrap = False
    tr = ccrs.AzimuthalEquidistant(central_longitude=cenlon,
                                   central_latitude=cenlat)
    edgefiles = [os.path.join(os.path.split(__file__)[0], 'edges_north.txt')]
    fbase = 'azeq_north'
    if notext:
        title = ''
    else:
        title = "NASA TESS's View\nof the Northern\nHemisphere"
    # turn off ecliptic coordinates since it doesn't matter
    ecliptic_coords = False
    # for now, don't try to put Kepler/K2 on the hemisphere maps
    addkepler = False
else:
    raise Exception(f'Unidentified hemisphere option: {hemisphere}')


# flag indicating we're just testing things
test = False
# create the output figure
makefig = True
# the output figure in high or "full" resolution
highres = True
fullres = False

# which color bar to use
color = 'gray'
if color == 'blue':
    cc = '_blue'
else:
    cc = ''

# save the output figure
savefig = True
# save every sector image for a gif in a subdirectory
makegif = False
if makegif:
    figdir = os.path.join(figdir, f'gif_{fbase}{cc}')
# use a transparent background instead of white
transparent = True
# the output figure file name
if transparent:
    if makegif:
        figdir += '_transp'
    fname = f'transp_{fbase}{cc}.png'
else:
    fname = f'{fbase}{cc}.png'
savefile = os.path.join(figdir, fname)

# get font sizes right for the output image size
if highres:
    fscl = 95
    if hemisphere == 'both':
        xinch = 150
        yinch = 75
    else:
        xinch = 100
        yinch = 100
    fsz = int(160 * fscl/100.)
    sfsz = int(175 * fscl/100.)
    tfsz = int(200 * fscl/100.)
elif fullres:
    fscl = 400
    if hemisphere == 'both':
        xinch = 600
        yinch = 300
    else:
        xinch = 400
        yinch = 400
    fsz = int(160 * fscl/100.)
    sfsz = int(175 * fscl/100.)
    tfsz = int(200 * fscl/100.)
else:
    if hemisphere == 'both':
        xinch = 12
        yinch = 6
    else:
        xinch = 8
        yinch = 8
    fsz = 12
    sfsz = 13
    tfsz = 15

# for creating the empirical corner models
xxs, yys, dats, ccds = [], [], [], []

# create the figure
if makefig:
    fig = plt.figure(figsize=(xinch, yinch))
    if hemisphere == 'both' and platecarree:
        ax = plt.axes([0.0, 0.0, 1.0, 1.0], projection=tr)
    else:
        # 1% border on all sides
        ax = plt.axes([0.01, 0.01, 0.98, 0.98], projection=tr)
        # ax = plt.axes([0.0, 0.0, 1.0, 1.0], projection=tr)
    if not test:
        # remove the outline circle of the globe
        ax.spines['geo'].set_linewidth(0)
    # set transparency
    if transparent:
        ax.background_patch.set_alpha(0)
    # the data coordinates are lat/lon in a grid
    data_tr = ccrs.PlateCarree()
    # load the edges of all the outer CCDs and invisibly plot them
    # so that after just 1 sector, the plot isn't artificially
    # zoomed to just that one sector.
    for edgefile in edgefiles:
        elat, elon = np.loadtxt(edgefile, unpack=True)
        if not test:
            plt.scatter(elon, elat, c='w', alpha=0.01, zorder=-5,
                        marker='.', s=1, transform=data_tr)
        else:
            pass
            # plt.scatter(elon, elat, c='r', alpha=1, zorder=5,
            #             s=20, transform=data_tr)

tois = pd.read_csv(toilist)

if useexoplots:
    coords = SkyCoord(ra=tois['ra'], dec=tois['dec'], unit=(u.deg, u.deg))
else:
    coords = SkyCoord(ra=tois['RA'], dec=tois['Dec'], unit=(u.hourangle, u.deg))
ecoords = coords.transform_to(BarycentricTrueEcliptic)
eclon = ecoords.lon.value * 1
eclat = ecoords.lat.value * 1
# transform to ecliptic coordinates if desired
if ecliptic_coords:
    lon = ecoords.lon.value * 1
    lat = ecoords.lat.value * 1
else:
    lon = coords.ra.value * 1
    lat = coords.dec.value * 1

# lon must be between -180 and 180 instead of 0 to 360
lon -= 180.
# because in astronomy images, plots have east on the left,
# switch east and west
lon *= -1.

if useexoplots:
    cands = (tois['disposition'] == 'Candidate') & (tois['facility'] == 'TESS')
    conf = (tois['disposition'] == 'Confirmed') & (tois['facility'] == 'TESS')
else:
    for ii in np.arange(tois.shape[0]):
        if type(tois.loc[ii, 'TFOPWG Disposition']) is not str:
            tois.loc[ii, 'TFOPWG Disposition'] = 'PC'

    cands = (tois['TFOPWG Disposition'] == 'APC') | \
            (tois['TFOPWG Disposition'] == 'PC')
    conf = (tois['TFOPWG Disposition'] == 'CP') | \
           (tois['TFOPWG Disposition'] == 'KP')

"""
norths = []
souths = []

for ii in np.arange(tois.shape[0]):
    secs = np.array(tois.loc[ii, 'Sectors'].split(',')).astype(int)
    north = False
    south = False
    for isec in secs:
        if isec <=13 or (isec >= 27 and isec <=39):
            south = True
        if (isec >=14 and isec <= 26) or (isec >=40 and isec <= 41) or 
                (isec >= 47):
            north = True
    norths.append(north)
    souths.append(south)

norths = np.array(norths)
souths = np.array(souths)
"""

norths = eclat > 0.73
# XXX: the gap between these is where the CCD gap was
souths = eclat < 0.37

if hemisphere == 'north':
    selcn = cands & norths
    selcf = conf & norths
elif hemisphere == 'south':
    selcn = cands & souths
    selcf = conf & souths
else:
    selcn = cands
    selcf = conf

if highres:
    cansz = 300
    consz = 1300
else:
    cansz = 5
    consz = 15

plt.scatter(lon[selcn], lat[selcn], c='#F1A93B', alpha=1, zorder=1,
            marker='o', s=cansz, transform=data_tr)

plt.scatter(lon[selcf], lat[selcf], c='#7BF9FC', alpha=1, zorder=1,
            marker='o', s=consz, transform=data_tr)

if savefig:
    plt.savefig(os.path.join(figdir, hemisphere+'.png'),
                transparent=transparent)
