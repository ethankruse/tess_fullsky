#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
from glob import glob
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from astropy.io import fits
from astropy.wcs import WCS
import cartopy.crs as ccrs
from truncate import truncate_colormap
from clean import clean_corner
from astropy.coordinates import SkyCoord, BarycentricTrueEcliptic

##################################################################
# Configuration parameters

datadir = os.path.join(os.path.split(__file__)[0], 'data')
figdir = os.path.join(os.path.split(__file__)[0], 'figs')
cornerdir = os.path.join(os.path.split(__file__)[0], 'corners')

# options are 'north', 'south', or 'both'
hemisphere = 'south'
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

# minimum and maximum flux for the colorbar
vmin = 150
vmax = 901.

# do we need to create the empirical corner glow correction for a sector?
makecorner = False
cornersec = 38

# remove the corner glow from the final image
remove_corner_glow = True
# make a plot of the corner glow for every CCD to check how removal is working
corner_glow_plot = False

# manual adjustments to the strength of corner glow corrections
adjfile = os.path.join(cornerdir, 'adjustments.txt')

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
makegif = True
if makegif:
    figdir = os.path.join(figdir, f'gif_{fbase}{cc}')
# use a transparent background instead of white
transparent = False
# the output figure file name
if transparent:
    if makegif:
        figdir += '_transp'
    fname = f'transp_{fbase}{cc}.png'
else:
    fname = f'{fbase}{cc}.png'
savefile = os.path.join(figdir, fname)

# credit text in lower left
if notext:
    credit = ''
else:
    credit = 'By Ethan Kruse\n@ethan_kruse'

# print the dates of the images in the lower right
printdate = True
if notext:
    printdate = False

# dates the sectors started and ended for when we print these
secstarts = {1: 'Jul 2018', 2: 'Aug 2018', 3: 'Sep 2018', 4: 'Oct 2018',
             5: 'Nov 2018', 6: 'Dec 2018', 7: 'Jan 2019', 8: 'Feb 2019',
             9: 'Feb 2019', 10: 'Mar 2019', 11: 'Apr 2019', 12: 'May 2019',
             13: 'Jun 2019', 14: 'Jul 2019', 15: 'Aug 2019', 16: 'Sep 2019', 
             17: 'Oct 2019', 18: 'Nov 2019', 19: 'Nov 2019', 20: 'Dec 2019',
             21: 'Jan 2020', 22: 'Feb 2020', 23: 'Mar 2020', 24: 'Apr 2020', 
             25: 'May 2020', 26: 'Jun 2020', 27: 'Jul 2020', 28: 'Jul 2020',
             29: 'Aug 2020', 30: 'Sep 2020', 31: 'Oct 2020', 32: 'Nov 2020',
             33: 'Dec 2020', 34: 'Jan 2021', 35: 'Feb 2021', 36: 'Mar 2021', 
             37: 'Apr 2021', 38: 'Apr 2021', 39: 'May 2021'}
secends = {1: 'Aug 2018', 2: 'Sep 2018', 3: 'Oct 2018', 4: 'Nov 2018',
           5: 'Dec 2018', 6: 'Jan 2019', 7: 'Feb 2019', 8: 'Feb 2019',
           9: 'Mar 2019', 10: 'Apr 2019', 11: 'May 2019', 12: 'Jun 2019',
           13: 'Jul 2019', 14: 'Aug 2019', 15: 'Sep 2019', 16: 'Oct 2019', 
           17: 'Nov 2019', 18: 'Nov 2019', 19: 'Dec 2019', 20: 'Jan 2020',
           21: 'Feb 2020', 22: 'Mar 2020', 23: 'Apr 2020', 24: 'May 2020', 
           25: 'Jun 2020', 26: 'Jul 2020', 27: 'Jul 2020', 28: 'Aug 2020',
           29: 'Sep 2020', 30: 'Oct 2020', 31: 'Nov 2020', 32: 'Dec 2020',
           33: 'Jan 2021', 34: 'Feb 2021', 35: 'Mar 2021', 36: 'Apr 2021', 
           37: 'Apr 2021', 38: 'May 2021', 39: 'Jun 2021'}

##################################################################

cnorm = colors.LogNorm(vmin=vmin, vmax=vmax)
if color == 'gray':
    # set up our custom colormap, which is a subset of the matplotlib map 'gray'.
    # we use truncate_colormap() to remove the blackest part of the map
    # so that even the dark areas show up against a pure black background.
    cmap = 'gray'
    # use only the latter part of the original colormap
    cmap = truncate_colormap(plt.get_cmap(cmap), minval=0.18, maxval=1.0)
elif color == 'blue':
    from matplotlib.colors import LinearSegmentedColormap
    carr = np.loadtxt('blue_cbar.txt')
    cmap = LinearSegmentedColormap.from_list('my_cmap', carr)
else:
    raise Exception('not recognized color')

if fullres and highres:
    raise Exception('Must choose either full or high resolution.')

# download the necessary files if they don't already exist
download = os.path.join(datadir, 'download.sh')
with open(download, 'r') as ff:
    for line in ff.readlines():
        fname = os.path.join(datadir, line.split()[5])
        if not os.path.exists(fname):
            subprocess.run(line, shell=True, cwd=datadir)
            
if addkepler:
    download2 = os.path.join(datadir, 'kepler_k2.sh')
    with open(download2, 'r') as ff:
        for line in ff.readlines():
            fname = os.path.join(datadir, line.split()[-2][1:-1])
            if not os.path.exists(fname):
                subprocess.run(line, shell=True, cwd=datadir)

if addkepler:
    allfiles = glob(os.path.join(datadir, '*fits'))
else:
    allfiles = glob(os.path.join(datadir, 'tess*fits'))
allfiles.sort()

# only grab the correct hemispheres of data
files = []
kepfiles = []
for ifile in allfiles:
    # grab the sector from the file name
    if os.path.split(ifile)[1][0] == 't':
        fsec = int(os.path.split(ifile)[1].split('-')[1][1:])
        # decide if we want to use it
        if ((fsec < 14) or (fsec > 26)) and hemisphere in ['both', 'south']:
            files.append(ifile)
        elif 13 < fsec < 27 and hemisphere in ['both', 'north']:
            files.append(ifile)
    else:
        if addkepler and hemisphere == 'both':
            kepfiles.append(ifile)
        else:
            raise Exception('Should not have non-TESS files in the list')

# make sure the output directory exists
if not os.path.exists(figdir) and makefig:
    os.makedirs(figdir, exist_ok=True)

# if we're creating the corner glow adjustments, only use that sector's files
# and make sure we're not trying to make a final output
if makecorner:
    files = glob(os.path.join(datadir, f'*s00{cornersec:02d}-*fits'))
    files.sort()
    makefig = False
    savefig = False
    makegif = False
    remove_corner_glow = True
    test = False
    kepfiles = []

# remove any previous images in the gif subdirectory
if makegif:
    prev = glob(os.path.join(figdir, '*png'))
    for iprev in prev:
        os.remove(iprev)

# anything we want to test
if test:
    files = glob(os.path.join(datadir, f'*s0037-3-3*fits'))
    #files = glob(os.path.join(datadir, f'*s0035-2-1*fits'))
    files += glob(os.path.join(datadir, f'*s0038-3-3*fits'))
    
    #files += glob(os.path.join(datadir, f'*s0011-1-3*fits'))
    #files += glob(os.path.join(datadir, f'*s0011-3-2*fits'))

    #files += glob(os.path.join(datadir, f'*s0007-2-1*fits'))
    #files += glob(os.path.join(datadir, f'*s0007-2-4*fits'))
    files.sort()
    
    
    kepfiles = []
    #kepfiles = glob(os.path.join(datadir, f'k*c08*fits'))
    #kepfiles += glob(os.path.join(datadir, f'k*c13*fits'))
    kepfiles.sort()


def grab_sector(sector, frac=0.95):
    """
    Print the download scripts for all the FFIs from a certain fraction of
    the way through a sector.

    Assumes that the tesscurl_sector_X_ffic.sh file already exists in the
    data directory.

    Parameters
    ----------
    sector : int
    frac : float
        What fraction (0-1) of the way through the sector do we want to get FFIs
        to use in this graphic

    Returns
    -------
    List[str]
    """
    # get time-ordered list of all FFI files in the sector
    sfile = os.path.join(datadir, f'tesscurl_sector_{sector}_ffic.sh')
    with open(sfile, 'r') as sff:
        lines = sff.readlines()
    lines.sort()

    # go to the right fraction of the way through the sector
    fthru = int(len(lines)*frac)
    ll = lines[fthru]
    # figure out what timestamp this file has
    date = ll.split('tess')[1].split('-')[0]
    # grab all files with that same timestamp
    ret = [sline for sline in lines if date in sline]
    # make sure we have one image per CCD and print/return them
    assert len(ret) == 16
    for iret in ret:
        print(iret[:-1])
    return ret


# load the manual corner glow adjustments
bsec, bcam, bccd, badj = np.loadtxt(adjfile, unpack=True, ndmin=2,
                                    delimiter=',', dtype=float)
bsec = bsec.astype(int)
bcam = bcam.astype(int)
bccd = bccd.astype(int)
# package these up
adjustments = (bsec, bcam, bccd, badj)

plt.close('all')

# get font sizes right for the output image size
if highres:
    fscl = 80
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
            #plt.scatter(elon, elat, c='r', alpha=1, zorder=5,
            #            s=20, transform=data_tr)

    # add the labels
    plt.text(0.02, 0.02, credit, transform=fig.transFigure,
             ha='left', va='bottom', multialignment='left',
             fontsize=fsz, fontname='Carlito')
    plt.text(0.02, 0.98, title, transform=fig.transFigure,
             ha='left', va='top', multialignment='left',
             fontsize=tfsz, fontname='Carlito')
    

def rebin(arr, new_shape):
    """Rebin 2D array arr to shape new_shape by averaging."""
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)


for ict, ifile in enumerate(kepfiles):
    now = datetime.now().strftime("%d.%m.%Y %H:%M:%S") 
    print(f'{now}. Processing Kepler/K2 image {ict+1} of {len(kepfiles)}.')
    with fits.open(ifile) as ffi:
        if 'campaign' in ffi[0].header:
            camp = ffi[0].header['campaign']
        else:
            camp = 'kep'
        # go through each CCD in a Kepler FFI
        for ii in np.arange(1, len(ffi)):
            wcs = WCS(ffi[ii].header)
            data = ffi[ii].data * 1
            # convert from max Kepler counts to max TESS counts
            data *= 0.6
            
            # get the coordinates of every data point in CCD coordinates
            xinds = np.arange(-0.5, data.shape[0]-0.4)
            yinds = np.arange(-0.5, data.shape[1]-0.4)
            mesh = np.meshgrid(xinds, yinds, indexing='ij')
            
            # transofrm to actual sky coordinates
            lon, lat = wcs.all_pix2world(mesh[1].flatten(), mesh[0].flatten(), 0)
            lon = lon.reshape(mesh[0].shape)
            lat = lat.reshape(mesh[1].shape)
            
            # detect bad modules
            if (lon[0,0] == 0.5) and (lat[0,0] == 0.5):
                continue
            
            # chop out unexposed rows/columns
            data = data[20:1044, 12:1112]
            lon = lon[20:1045, 12:1113]
            lat = lat[20:1045, 12:1113]
        
            # transform to ecliptic coordinates if desired
            if ecliptic_coords:
                icrs = SkyCoord(ra=lon, dec=lat, frame='icrs', unit='deg')
                ecliptic = icrs.transform_to(BarycentricTrueEcliptic)
                lon = ecliptic.lon.value * 1
                lat = ecliptic.lat.value * 1
        
            # lon must be between -180 and 180 instead of 0 to 360
            lon -= 180.
            # because in astronomy images, plots have east on the left,
            # switch east and west
            lon *= -1.
        
            # rebin from Kepler's 4 arcsec to 16 arcsec pixels, closer to TESS
            data = rebin(data, (data.shape[0]//4, data.shape[1]//4))
            lat = lat[::4, ::4]
            lon = lon[::4, ::4]
            
            # make adjustments to match the TESS intersections and each other
            # campaigns 10 and 11 were split into two pieces and these ints
            # work (111/112 and 101/102)
            if camp in [0, 6, 15, 19]:
                data -= 10
            elif camp in [2, 4, 5, 8, 13, 18]:
                data -= 20
            elif camp in [17]:
                data -= 30
            elif camp in [111, 16]:
                data -= 40
            elif camp in [3, 12]:
                data -= 50.
            
            if makefig:
                # for wraparounds:
                if wrap and lon.max() > cenlon + 178 and lon.min() < cenlon - 178:                    
                    # find the problem areas that wrap around in longitude
                    bad = ((np.abs(lon[:-1,:-1] - lon[:-1,1:]) > 355.)|
                           (np.abs(lon[:-1,:-1] - lon[1:,:-1]) > 355.)|
                           (np.abs(lon[:-1,:-1] - lon[1:,1:]) > 355.))
                    # mask them and just don't plot these pixels
                    maskeddata = np.ma.masked_where(bad, data)
                    plt.pcolormesh(lon, lat, maskeddata, norm=cnorm, alpha=1,
                                   transform=data_tr, cmap=cmap)
                else:
                    # plot the actual image from this CCD
                    plt.pcolormesh(lon, lat, data, norm=cnorm, alpha=1,
                                   transform=data_tr, cmap=cmap)
                    
# save the plot after each sector for the gif
if makegif and savefig and makefig and len(kepfiles) > 0:
    if transparent:
        outfig = os.path.join(figdir, f'transp_img{0:04d}.png')
    else:
        outfig = os.path.join(figdir, f'img{0:04d}.png')
    plt.savefig(outfig, transparent=transparent)

# loop through every image and create the mosaic
for ii, ifile in enumerate(files):
    with fits.open(ifile) as ff:
        wcs = WCS(ff[1].header)
        data = ff[1].data * 1

        # get the coordinates of every data point in CCD coordinates
        xinds = np.arange(-0.5, data.shape[0]-0.4)
        yinds = np.arange(-0.5, data.shape[1]-0.4)
        mesh = np.meshgrid(xinds, yinds, indexing='ij')

        # transofrm to actual sky coordinates
        lon, lat = wcs.all_pix2world(mesh[1].flatten(), mesh[0].flatten(), 0)
        lon = lon.reshape(mesh[0].shape)
        lat = lat.reshape(mesh[1].shape)

        # chop out unexposed rows/columns
        # this should always be true, leaving a 2048x2048 image
        data = data[:2048, 44:2092]
        lon = lon[:2049, 44:2093]
        lat = lat[:2049, 44:2093]
        
        # transform to ecliptic coordinates if desired
        if ecliptic_coords:
            icrs = SkyCoord(ra=lon, dec=lat, frame='icrs', unit='deg')
            ecliptic = icrs.transform_to(BarycentricTrueEcliptic)
            lon = ecliptic.lon.value * 1
            lat = ecliptic.lat.value * 1

        # lon must be between -180 and 180 instead of 0 to 360
        lon -= 180.
        # because in astronomy images, plots have east on the left,
        # switch east and west
        lon *= -1.

        # what image is this
        icam = ff[1].header['camera']
        iccd = ff[1].header['ccd']
        isec = int(ifile.split('-s0')[1][:3])
        now = datetime.now().strftime("%d.%m.%Y %H:%M:%S") 
        print(f'{now}. Processing image {ii+1} of {len(files)}: Sector {isec} '
              f'Camera {icam} CCD {iccd}')

        if remove_corner_glow:
            # create the empirical corner glow model
            if makecorner:
                xs, ys, dat = clean_corner(data, cornerdir, adjustments,
                                           cleanplot=corner_glow_plot,
                                           create=True, ccd=iccd, sec=isec,
                                           cam=icam)
                # save the coordinates and smoothed corner
                xxs.append(xs)
                yys.append(ys)
                dats.append(dat)
                ccds.append(iccd)
            # remove the corner glow from this CCD
            else:
                data = clean_corner(data, cornerdir, adjustments,
                                    cleanplot=corner_glow_plot, ccd=iccd,
                                    sec=isec, cam=icam)
                

        # some special processing to avoid problem areas.
        # these are all in overlap zones, so removing contaminated chunks
        # of these images doesn't affect the overall footprint and makes for a
        # cleaner final image.
        if isec == 12 and icam == 4 and iccd == 3:
            # avoid the contamination glow from the slightly off-camera bright
            # star Canopus
            data = data[300:, :]
            lat = lat[300:, :]
            lon = lon[300:, :]
        elif isec == 11 and icam == 3 and iccd == 3:
            # a reflection in the CVZ only in S11
            data = data[450:, :]
            lat = lat[450:, :]
            lon = lon[450:, :]
        elif isec == 9 and icam == 3 and iccd == 3:
            # a reflection in the CVZ only in S9
            data = data[400:, :]
            lat = lat[400:, :]
            lon = lon[400:, :]
        elif isec == 19 and icam == 3 and iccd == 3:
            # a reflection in the CVZ only in S19
            data = data[360:, :]
            lat = lat[360:, :]
            lon = lon[360:, :]
        elif isec == 21 and icam == 2 and iccd == 1:
            # a reflection in the CVZ only in S21
            data = data[225:, :]
            lat = lat[225:, :]
            lon = lon[225:, :]
        elif isec == 23 and icam == 3 and iccd == 4:
            # bad corner glow
            data[:200, -100:] = np.nan
        elif isec == 24 and icam == 2 and iccd == 2:
            data[:330, 1850:] = np.nan
        elif isec == 24 and icam == 3 and iccd == 3:
            data[:200, :700] = np.nan
        elif isec == 24 and icam == 4 and iccd == 3:
            data[:150, 1700:] = np.nan
            data[:500, :300] = np.nan
        elif isec in [24, 30] and icam == 4 and iccd == 4:
            data[:400, 1600:] = np.nan
        elif isec == 26 and icam == 1 and iccd == 1:
            data[:350, :350] = np.nan
        elif isec == 26 and icam == 3 and iccd == 3:
            data[:200, :500] = np.nan
        elif isec in [26, 28, 29] and icam == 4 and iccd == 3:
            data[:700, :500] = np.nan
        elif isec == 27 and icam == 1 and iccd == 4:
            data[:700, -600:] = np.nan
        elif isec == 27 and icam == 1 and iccd == 2:
            data[:120, -120:] = np.nan
        elif isec == 27 and icam == 3 and iccd == 3:
            data[:150, 300:550] = np.nan
        elif isec == 29 and icam == 1 and iccd == 1:
            data[:200, :350] = np.nan
        elif isec == 30 and icam == 1 and iccd == 4:
            data[350:650, 1900:] = np.nan
        elif isec == 5 and icam == 2 and iccd == 1:
            data[:200, :350] = np.nan
        elif isec == 28 and icam == 2 and iccd == 1:
            data[450:650, :100] = np.nan
        elif isec == 28 and icam == 3 and iccd == 2:
            data[250:450, 1900:] = np.nan
        elif isec == 13 and icam == 4 and iccd == 3:
            data[:400, :125] = np.nan
        elif isec == 13 and icam == 2 and iccd == 1:
            data[:100, 500:665] = np.nan
        elif isec == 11 and icam == 2 and iccd == 1:
            data[:150, 300:665] = np.nan        
        elif isec == 10 and icam == 2 and iccd == 1:
            data[:250, 100:350] = np.nan 
        elif isec == 9 and icam == 2 and iccd == 1:
            data[:100, :100] = np.nan 
        elif isec == 9 and icam == 3 and iccd == 2:
            data[300:700, 1750:] = np.nan 
        elif isec == 8 and icam == 3 and iccd == 2:
            data[:100, 1700:1950] = np.nan 
        elif isec == 31 and icam == 1 and iccd == 1:
            data[50:300, 150:400] = np.nan
        elif isec == 31 and icam == 1 and iccd == 2:
            data[50:300, 1850:] = np.nan
        elif isec == 31 and icam == 2 and iccd == 3:
            data[150:300, :100] = np.nan
        elif isec == 31 and icam == 3 and iccd == 4:
            data[:300, 1300:1700] = np.nan
        elif isec == 31 and icam == 4 and iccd == 4:
            data[:500, 1400:] = np.nan
        elif isec == 32 and icam == 1 and iccd == 4:
            data[500:700, 1800:] = np.nan
            data[:300, 1200:1700] = np.nan
        elif isec == 32 and icam == 2 and iccd == 3:
            data[:150, 200:700] = np.nan
        elif isec == 32 and icam == 2 and iccd == 4:
            # remove saturated columns
            data[:2010, 772] = data[:2010, 771]
            data[:2010, 773] = data[:2010, 774]
        elif isec == 32 and icam == 3 and iccd == 4:
            data[:150, 1200:1500] = np.nan
            data[:200, 1900:] = np.nan
        elif isec == 33 and icam == 1 and iccd == 1:
            data[:150, 450:750] = np.nan
        elif isec == 33 and icam == 2 and iccd == 4:
            data[:600, 1750:] = np.nan
        elif isec == 33 and icam == 1 and iccd == 3:
            data[:350, :300] = np.nan
        elif isec == 33 and icam == 1 and iccd == 2:
            data[:200, 1800:] = np.nan
            data[:140, 1440:1550] = np.nan
        elif isec == 33 and icam == 3 and iccd == 1:
            data[:40, :60] = np.nan
        elif isec == 17 and icam == 3 and iccd == 4:
            data[:200, 1750:] = np.nan
        elif isec == 20 and icam == 3 and iccd == 3:
            data[200:450, :200] = np.nan
        elif isec == 20 and icam == 4 and iccd == 4:
            data[:200, 1700:] = np.nan
        elif isec == 23 and icam == 4 and iccd == 3:
            data[:150, :300] = np.nan
        elif isec == 26 and icam == 1 and iccd == 3:
            data[20:130, 250:400] = np.nan
            data[300:450, 400:650] = np.nan
        elif isec == 34 and icam == 1 and iccd == 4:
            data[:600, 1400:] = np.nan
        elif isec == 34 and icam == 1 and iccd == 1:
            data[:600, :700] = np.nan
        elif isec == 34 and icam == 2 and iccd == 4:
            data[:700, 1500:] = np.nan
        elif isec == 35 and icam == 2 and iccd == 1:
            data[:600, :700] = np.nan
            data[:1900, :300] = np.nan
        elif isec == 35 and icam == 3 and iccd == 2:
            data[:400, 1300:] = np.nan
            data[:800, 1900:] = np.nan
            data[:, 1900:] = np.nan
        elif isec == 35 and icam == 3 and iccd == 3:
            data[:, :250] = np.nan
        elif isec == 36 and icam == 2 and iccd == 1:
            data[:300, :600] = np.nan
        elif isec == 36 and icam == 3 and iccd == 2:
            data[:800, 1500:] = np.nan
        elif isec == 36 and icam == 3 and iccd == 3:
            data[:400, :400] = np.nan
        elif isec == 37 and icam == 1 and iccd == 1:
            data[:250, :300] = np.nan
        elif isec == 37 and icam == 2 and iccd == 1:
            data[:275, :900] = np.nan
        elif isec == 37 and icam == 3 and iccd == 2:
            data[:300, 1700:] = np.nan
        elif isec == 37 and icam == 3 and iccd == 3:
            data[250:350, 400:600] = np.nan
        elif isec == 37 and icam == 2 and iccd == 2:
            data[:600, 1400:] = np.nan
        elif isec == 37 and icam == 3 and iccd == 1:
            data[:200, :500] = np.nan
        elif isec == 37 and icam == 2 and iccd == 4:
            data[250:360, 1500:1625] = np.nan
        elif isec == 38 and icam == 2 and iccd == 4:
            data[:450, 1700:] = np.nan
        elif isec == 38 and icam == 3 and iccd == 3:
            data[:250, :300] = np.nan
        elif isec == 38 and icam == 4 and iccd == 2:
            data[:400, 1800:] = np.nan
        elif isec == 38 and icam == 2 and iccd == 1:
            data[:600, :700] = np.nan
        elif isec == 38 and icam == 3 and iccd == 2:
            data[:500, 1500:] = np.nan
            
        # remove weird saturated columns that don't have obvious sources
        if isec == 26 and icam == 3 and iccd == 3:
            data[:, 1195] = data[:, 1194]
            data[:, 1196] = data[:, 1197]
        if isec == 26 and icam == 2 and iccd == 4:
            data[:, 507] = data[:, 506]
            data[:, 508] = data[:, 509]
        
        # this camera in this sector is too bright and doesn't match the rest
        if isec == 24 and icam == 4:
            if iccd in [3, 4]:
                data -= 70
            else:
                data -= 35
        if isec == 25 and icam == 4:
            if iccd in [3, 4]:
                data -= 50
            else:
                data -= 25
        if isec == 26 and icam == 4 and iccd in [3, 4]:
            data -= 35
        if isec == 15 and icam == 1 and iccd in [2, 3]:
            data += 15
        if isec == 16 and icam == 1 and iccd in [2]:
            data += 20
        if isec == 27 and icam == 1 and iccd in [3, 4]:
            data += 23
        if isec == 28 and icam == 1 and iccd in [3, 4]:
            data += 15
        if isec in [27, 28] and icam == 1 and iccd in [1, 2]:
            data += 15
        if isec == 1 and icam == 2 and iccd == 2:
            data -= 30
        if isec == 6 and icam == 1 and iccd == 1:
            data -= 15
        if isec == 7 and icam == 1 and iccd == 1:
            data -= 15
        if isec == 7 and icam == 1 and iccd == 4:
            data -= 10
        if isec == 38 and icam == 2 and iccd == 1:
            data += 15
        if isec == 38 and icam == 3 and iccd == 3:
            data += 10
        """
        if isec == 35 and icam == 2 and iccd == 1:
            data += 15
        if isec == 35 and icam == 3 and iccd == 2:
            data += 15
        if isec == 35 and icam == 3 and iccd == 3:
            data += 15
        """
        if corner_glow_plot:
            plt.figure()
            plt.imshow(data, norm=cnorm, cmap=cmap)

        if makefig:
            if ii == 0:
                sectxt = f'Sector {isec}\n{secstarts[isec]}\u2013{secends[isec]}'
                if printdate:
                    text = plt.text(0.98, 0.02, sectxt,
                                    transform=fig.transFigure, ha='right',
                                    va='bottom', multialignment='right',
                                    fontsize=sfsz, fontname='Carlito')
                ssec = isec
            # for wraparounds:
            if wrap and lon.max() > cenlon + 178 and lon.min() < cenlon - 178:                    
                # find the problem areas that wrap around in longitude
                bad = ((np.abs(lon[:-1,:-1] - lon[:-1,1:]) > 355.)|
                       (np.abs(lon[:-1,:-1] - lon[1:,:-1]) > 355.)|
                       (np.abs(lon[:-1,:-1] - lon[1:,1:]) > 355.))
                # mask them and just don't plot these pixels
                maskeddata = np.ma.masked_where(bad, data)
                plt.pcolormesh(lon, lat, maskeddata, norm=cnorm, alpha=1,
                               transform=data_tr, cmap=cmap)

            else:
                # plot the actual image from this CCD
                plt.pcolormesh(lon, lat, data, norm=cnorm, alpha=1,
                               transform=data_tr, cmap=cmap)

        # if we're starting to plot a new sector, update the date
        if (ii % 16) == 0 and ii > 0 and printdate and makefig:
            text.remove()
            if hemisphere == 'both' or isec < 27:
                sectxt = f'Sectors {ssec}\u2013{isec}\n{secstarts[ssec]}\u2013{secends[isec]}'
            elif hemisphere == 'south':
                if isec == 27:
                    sectxt = f'Sectors 1\u201313; 27\n{secstarts[1]}\u2013{secends[13]}\n{secstarts[27]}\u2013{secends[27]}'
                else:
                    sectxt = f'Sectors 1\u201313; 27\u2013{isec}\n{secstarts[1]}\u2013{secends[13]}\n{secstarts[27]}\u2013{secends[isec]}'
            else:
                raise Exception('need to handle this')
            text = plt.text(0.98, 0.02, sectxt, transform=fig.transFigure,
                            ha='right', va='bottom', multialignment='right',
                            fontsize=sfsz, fontname='Carlito')
        # save the plot after each sector for the gif
        if makegif and savefig and makefig and ii > 0 and ((ii+1) % 16) == 0:
            if transparent:
                outfig = os.path.join(figdir, f'transp_img{(ii+1)//16:04d}.png')
            else:
                outfig = os.path.join(figdir, f'img{(ii+1)//16:04d}.png')
            plt.savefig(outfig, transparent=transparent)

# save the figure if we haven't already as part of the gif
# don't overwrite any previous file, increment the number
if makefig and savefig and not makegif:
    inum = 1
    orig = savefile
    while os.path.exists(savefile):
        savefile = (os.path.splitext(orig)[0] + f'{inum}' +
                    os.path.splitext(orig)[1])
        inum += 1
    plt.savefig(savefile, transparent=transparent)
    
now = datetime.now().strftime("%d.%m.%Y %H:%M:%S") 
print(f'{now}. Finished.')

# if we're creating the empirical corner model, look at all the corners from
# this sector together
if makecorner:
    if not os.path.exists(cornerdir):
        os.makedirs(cornerdir, exist_ok=True)

    plt.close('all')
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    xs = xxs[0]
    ys = yys[0]

    stack = []
    for ii in np.arange(len(dats)):
        dat = dats[ii] * 1
        mccd = ccds[ii]

        if mccd == 2 or mccd == 4:
            corner = 3
        elif mccd == 1 or mccd == 3:
            corner = 1
        else:
            raise Exception('Bad CCD')

        # get the orientation of all the corner glow regions the same
        if corner == 3:
            dat = dat[:, ::-1]

        dx, dy = dat.shape
        # look at the 'baseline' noise level far enough from the corner to not
        # be strongly affected by the glow
        imin = np.median(dat[3*dx//4:, 3*dy//4:])
        # plot and save the baseline adjusted corner glow
        ax.plot_surface(xs, ys, dat - imin, cmap='gray',
                        linewidth=0, antialiased=False, alpha=0.2)
        stack.append(dat - imin)

    # get the empirical median corner glow shape for this sector and plot it in
    # color to compare against the stack of individual corners
    avg = np.median(np.dstack(stack), axis=-1)
    ax.plot_surface(xs, ys, avg, cmap='viridis',
                    linewidth=0, antialiased=False, alpha=0.4)
    plt.title(f'Sec {cornersec} Corners')

    # save the diagnostic plot and the median corner glow to be used in future
    # corrections
    ctxt = os.path.join(cornerdir, f'sector{cornersec:02d}.corner.txt')
    cfig = os.path.join(cornerdir, f'sector{cornersec:02d}.corner.png')
    np.savetxt(ctxt, avg)
    plt.savefig(cfig)


