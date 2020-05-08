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

# parameters that change depending on the hemisphere
if hemisphere == 'both':
    # coordinates at the center of the projection
    cenlon = 0.
    cenlat = 0.
    # if the projection has a dividing line where we have to wrap data around
    wrap = True
    # set up our desired map projection
    tr = ccrs.Mollweide(central_longitude=cenlon)
    # the coordinates of the corners of the CCDs
    edgefiles = [os.path.join(os.path.split(__file__)[0], 'edges_south.txt'), 
                 os.path.join(os.path.split(__file__)[0], 'edges_north.txt')]
    # what the output file name base should be
    fbase = 'mollweide'
    if ecliptic_coords:
        fbase += '_ecliptic'
    #  title text in upper left
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
    title = "NASA TESS's View\nof the Southern\nHemisphere"
    # turn off ecliptic coordinates since it doesn't matter
    ecliptic_coords = False
elif hemisphere == 'north':
    cenlon = -90.
    cenlat = 66.560708333333
    wrap = False
    tr = ccrs.AzimuthalEquidistant(central_longitude=cenlon,
                                   central_latitude=cenlat)
    edgefiles = [os.path.join(os.path.split(__file__)[0], 'edges_north.txt')]
    fbase = 'azeq_north'
    title = "NASA TESS's View\nof the Northern\nHemisphere"
    # turn off ecliptic coordinates since it doesn't matter
    ecliptic_coords = False
else:
    raise Exception(f'Unidentified hemisphere option: {hemisphere}')

# minimum and maximum flux for the colorbar
vmin = 150
vmax = 901.

# set up our custom colormap, which is a subset of the matplotlib map 'gray'.
# we use truncate_colormap() to remove the blackest part of the map
# so that even the dark areas show up against a pure black background.
cnorm = colors.LogNorm(vmin=vmin, vmax=vmax)
cmap = 'gray'
# use only the latter part of the original colormap
cmap = truncate_colormap(plt.get_cmap(cmap), minval=0.18, maxval=1.0)

# do we need to create the empirical corner glow correction for a sector?
makecorner = False
cornersec = 23

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
# the output figure in full resolution
highres = True
# save the output figure
savefig = True
# save every sector image for a gif in a subdirectory
makegif = True
if makegif:
    figdir = os.path.join(figdir, f'gif_{fbase}')
# use a transparent background instead of white
transparent = False
# the output figure file name
if transparent:
    fname = f'transp_{fbase}.png'
else:
    fname = f'{fbase}.png'
savefile = os.path.join(figdir, fname)

# credit text in lower left
credit = 'By Ethan Kruse\n@ethan_kruse'

# print the dates of the images in the lower right
printdate = True

# dates the sectors started and ended for when we print these
secstarts = {1: 'Jul 2018', 2: 'Aug 2018', 3: 'Sep 2018', 4: 'Oct 2018',
             5: 'Nov 2018', 6: 'Dec 2018', 7: 'Jan 2019', 8: 'Feb 2019',
             9: 'Feb 2019', 10: 'Mar 2019', 11: 'Apr 2019', 12: 'May 2019',
             13: 'Jun 2019', 14: 'Jul 2019', 15: 'Aug 2019', 16: 'Sep 2019', 
             17: 'Oct 2019', 18: 'Nov 2019', 19: 'Nov 2019', 20: 'Dec 2019',
             21: 'Jan 2020', 22: 'Feb 2020', 23: 'Mar 2020', 24: 'Apr 2020', 
             25: 'May 2020', 26: 'Jun 2020'}
secends = {1: 'Aug 2018', 2: 'Sep 2018', 3: 'Oct 2018', 4: 'Nov 2018',
           5: 'Dec 2018', 6: 'Jan 2019', 7: 'Feb 2019', 8: 'Feb 2019',
           9: 'Mar 2019', 10: 'Apr 2019', 11: 'May 2019', 12: 'Jun 2019',
           13: 'Jul 2019', 14: 'Aug 2019', 15: 'Sep 2019', 16: 'Oct 2019', 
           17: 'Nov 2019', 18: 'Nov 2019', 19: 'Dec 2019', 20: 'Jan 2020',
           21: 'Feb 2020', 22: 'Mar 2020', 23: 'Apr 2020', 24: 'May 2020', 
           25: 'Jun 2020', 26: 'Jul 2020'}

##################################################################

# download the necessary files if they don't already exist
download = os.path.join(datadir, 'download.sh')
with open(download, 'r') as ff:
    for line in ff.readlines():
        fname = os.path.join(datadir, line.split()[5])
        if not os.path.exists(fname):
            subprocess.run(line, shell=True, cwd=datadir)

allfiles = glob(os.path.join(datadir, '*fits'))
allfiles.sort()

# only grab the correct hemispheres of data
files = []
for ifile in allfiles:
    # grab the sector from the file name
    fsec = int(os.path.split(ifile)[1].split('-')[1][1:])
    # decide if we want to use it
    if fsec < 14 and hemisphere != 'north':
        files.append(ifile)
    elif fsec > 13 and hemisphere != 'south':
        files.append(ifile)

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

# remove any previous images in the gif subdirectory
if makegif:
    prev = glob(os.path.join(figdir, '*png'))
    for iprev in prev:
        os.remove(iprev)

# anything we want to test
if test:
    # files = files[187:188]
    #files = glob(os.path.join(datadir, f'*s0022-*fits'))
    #files += glob(os.path.join(datadir, f'*s0023*fits'))
    
    files = glob(os.path.join(datadir, f'*s0023-3-4*fits'))
    files += glob(os.path.join(datadir, f'*s0022-3-4*fits'))
    #files += glob(os.path.join(datadir, f'*s0018-2-4*fits'))
    files.sort()
    


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
        print(iret)
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
    fscl = 200
    if hemisphere == 'both':
        xinch = 300
        yinch = 150
    else:
        xinch = 200
        yinch = 200
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
                if corner_glow_plot:
                    plt.figure()
                    plt.imshow(data, norm=cnorm, cmap=cmap)

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

        if makefig:
            # create the figure. if testing, each CCD gets its own figure
            if ii == 0:
                fig = plt.figure(figsize=(xinch, yinch))
                # 1% border on all sides
                ax = plt.axes([0.01, 0.01, 0.98, 0.98], projection=tr)
                if not test:
                    # remove the outline circle of the globe
                    ax.outline_patch.set_linewidth(0)
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
                    plt.scatter(elon, elat, c='w', alpha=0.01, zorder=-5,
                                marker='.', s=1, transform=data_tr)
                    if test:
                        plt.scatter(elon, elat, c='r', alpha=1, zorder=5,
                                    s=20, transform=data_tr)

                # add the labels
                plt.text(0.02, 0.02, credit, transform=fig.transFigure,
                         ha='left', va='bottom', multialignment='left',
                         fontsize=fsz, fontname='Carlito')
                plt.text(0.02, 0.98, title, transform=fig.transFigure,
                         ha='left', va='top', multialignment='left',
                         fontsize=tfsz, fontname='Carlito')
                sectxt = f'Sector {isec}\n{secstarts[isec]}-{secends[isec]}'
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
        if (ii % 16) == 0 and ii > 0 and printdate:
            text.remove()
            sectxt = f'Sectors {ssec}-{isec}\n{secstarts[ssec]}-{secends[isec]}'
            text = plt.text(0.98, 0.02, sectxt, transform=fig.transFigure,
                            ha='right', va='bottom', multialignment='right',
                            fontsize=sfsz, fontname='Carlito')
        # save the plot after each sector for the gif
        if makegif and savefig and ii > 0 and ((ii+1) % 16) == 0:
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


