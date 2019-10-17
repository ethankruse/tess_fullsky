#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from astropy.io import fits
from astropy.wcs import WCS
import cartopy.crs as ccrs
from truncate import truncate_colormap
from clean import clean_corner

##################################################################
# Configuration parameters

datadir = os.path.join(os.path.split(__file__)[0], 'data')
figdir = os.path.join(os.path.split(__file__)[0], 'figs')
cornerdir = os.path.join(os.path.split(__file__)[0], 'corners')

# south ecliptic pole coordinates are 90, -66.560708333333
# coordinates at the center of the projection
cenlon = 90.
cenlat = -66.560708333333

# if the projection has a dividing line where we have to wrap data around
wrap = False

# set up our desired map projection
# tr = ccrs.Orthographic(central_longitude=cenlon,
#                         central_latitude=cenlat)
# tr = ccrs.Stereographic(central_longitude=cenlon,
#                         central_latitude=cenlat)
tr = ccrs.AzimuthalEquidistant(central_longitude=cenlon,
                               central_latitude=cenlat)

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
cornersec = 13

# remove the corner glow from the final image
remove_corner_glow = True
# make a plot of the corner glow for every CCD to check how removal is working
corner_glow_plot = False

# manual adjustments to the strength of corner glow corrections
adjfile = os.path.join(cornerdir, 'adjustments.txt')
# the coordinates of the corners of the CCDs
edgefile = os.path.join(os.path.split(__file__)[0], 'edges.txt')

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
    figdir = os.path.join(figdir, 'gif_azeq_label')
# use a transparent background instead of white
transparent = False
# the output figure file name
if transparent:
    fname = 'transp_ortho.png'
else:
    fname = 'ortho.png'
savefile = os.path.join(figdir, fname)

# credit text in lower left and title text in upper left
credit = 'By Ethan Kruse\n@ethan_kruse'
title = "NASA TESS's View\nof the Southern\nHemisphere"
# credit = ''
# title = ''

# print the dates of the images in the lower right
printdate = True

# dates the sectors started and ended for when we print these
secstarts = {1: 'Jul 2018', 2: 'Aug 2018', 3: 'Sep 2018', 4: 'Oct 2018',
             5: 'Nov 2018', 6: 'Dec 2018', 7: 'Jan 2019', 8: 'Feb 2019',
             9: 'Feb 2019', 10: 'Mar 2019', 11: 'Apr 2019', 12: 'May 2019',
             13: 'Jun 2019'}
secends = {1: 'Aug 2018', 2: 'Sep 2018', 3: 'Oct 2018', 4: 'Nov 2018',
           5: 'Dec 2018', 6: 'Jan 2019', 7: 'Feb 2019', 8: 'Feb 2019',
           9: 'Mar 2019', 10: 'Apr 2019', 11: 'May 2019', 12: 'Jun 2019',
           13: 'Jul 2019'}

##################################################################

# download the necessary files if they don't already exist
download = os.path.join(datadir, 'download.sh')
with open(download, 'r') as ff:
    for line in ff.readlines():
        fname = os.path.join(datadir, line.split()[5])
        if not os.path.exists(fname):
            subprocess.run(line, shell=True, cwd=datadir)

files = glob(os.path.join(datadir, '*fits'))
files.sort()

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
    pass


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
    inches = 200
    fsz = int(160 * inches/100.)
    sfsz = int(175 * inches/100.)
    tfsz = int(200 * inches/100.)
else:
    inches = 8
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

        # lon must be between -180 and 180 instead of 0 to 360
        lon -= 180.
        # because in astronomy images, plots have east on the left,
        # switch east and west
        lon *= -1.

        # what image is this
        icam = ff[1].header['camera']
        iccd = ff[1].header['ccd']
        isec = int(ifile.split('-s0')[1][:3])
        print(f'Processing image {ii+1} of {len(files)}: Sector {isec} Camera '
              f'{icam} CCD {iccd}')

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

        if makefig:
            # create the figure. if testing, each CCD gets its own figure
            if ii == 0 or test:
                if highres:
                    fig = plt.figure(figsize=(inches, inches))
                else:
                    fig = plt.figure()
                # 1% border on all sides
                ax = plt.axes([0.01, 0.01, 0.98, 0.98], projection=tr)
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
                elat, elon = np.loadtxt(edgefile, unpack=True)
                plt.scatter(elon, elat, c='w', alpha=0.01, zorder=-5,
                            marker='.', s=1, transform=data_tr)

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
            if wrap and lon.max() > cenlon + 120 and lon.min() < cenlon - 120:
                left = np.where(lon > cenlon + 120)
                lonleft = lon * 1
                lonleft[left] = cenlon - 180.
                plt.pcolormesh(lonleft, lat, data, norm=cnorm, alpha=1,
                               transform=data_tr, cmap=cmap)

                right = np.where(lon < cenlon - 120)
                lonright = lon * 1
                lonright[right] = cenlon + 180.
                plt.pcolormesh(lonright, lat, data, norm=cnorm, alpha=1,
                               transform=data_tr, cmap=cmap)
            # plot the actual image from this CCD
            else:
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
