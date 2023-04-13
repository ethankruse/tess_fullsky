import os
import subprocess
import warnings
from datetime import datetime
from glob import glob

import cartopy.crs as ccrs
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import BarycentricTrueEcliptic, Galactic, SkyCoord
from astropy.io import fits
from astropy.wcs import FITSFixedWarning, WCS
from matplotlib.image import imread

from clean import clean_corner
from truncate import truncate_colormap

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
galactic_coords = False
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
    if galactic_coords:
        cenlon += 180
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
    if galactic_coords:
        fbase += '_galactic'
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
    galactic_coords = False
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
    galactic_coords = False
    # for now, don't try to put Kepler/K2 on the hemisphere maps
    addkepler = False
else:
    raise Exception(f'Unidentified hemisphere option: {hemisphere}')

# minimum and maximum flux for the colorbar
vmin = 150
vmax = 901.

# do we need to create the empirical corner glow correction for a sector?
makecorner = False
cornersec = 62

# remove the corner glow from the final image
remove_corner_glow = True
# make a plot of the corner glow for every CCD to check how removal is working
corner_glow_plot = False

# manual adjustments to the strength of corner glow corrections
adjfile = os.path.join(cornerdir, 'adjustments.txt')

# save data/processing time by binning the images 2x2
binning = True

# flag indicating we're just testing things
test = False
# create the output figure
makefig = True
# the output figure in high or "full" resolution
highres = True
fullres = False

# which color bar to use
color = 'blue'
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
             37: 'Apr 2021', 38: 'Apr 2021', 39: 'May 2021', 40: 'Jun 2021',
             41: 'Jul 2021', 42: 'Aug 2021', 43: 'Sep 2021', 44: 'Oct 2021',
             45: 'Nov 2021', 46: 'Dec 2021', 47: 'Dec 2021', 48: 'Jan 2022',
             49: 'Feb 2022', 50: 'Mar 2022', 51: 'Apr 2022', 52: 'May 2022',
             53: 'Jun 2022', 54: 'Jul 2022', 55: 'Aug 2022', 56: 'Sep 2022',
             57: 'Sep 2022', 58: 'Oct 2022', 59: 'Nov 2022', 60: 'Dec 2022',
             61: 'Jan 2023', 62: 'Feb 2023', 63: 'Mar 2023', 64: 'Apr 2023',
             65: 'May 2023', 66: 'Jun 2023', 67: 'Jul 2023', 68: 'Jul 2023',
             69: 'Aug 2023', 70: 'Sep 2023', 71: 'Oct 2023', 72: 'Nov 2023',
             73: 'Dec 2023', 74: 'Jan 2024', 75: 'Jan 2024', 76: 'Feb 2024',
             77: 'Mar 2024', 78: 'Apr 2024', 79: 'May 2024', 80: 'Jun 2024',
             81: 'Jul 2024', 82: 'Aug 2024', 83: 'Sep 2024'}
secends = {1: 'Aug 2018', 2: 'Sep 2018', 3: 'Oct 2018', 4: 'Nov 2018',
           5: 'Dec 2018', 6: 'Jan 2019', 7: 'Feb 2019', 8: 'Feb 2019',
           9: 'Mar 2019', 10: 'Apr 2019', 11: 'May 2019', 12: 'Jun 2019',
           13: 'Jul 2019', 14: 'Aug 2019', 15: 'Sep 2019', 16: 'Oct 2019',
           17: 'Nov 2019', 18: 'Nov 2019', 19: 'Dec 2019', 20: 'Jan 2020',
           21: 'Feb 2020', 22: 'Mar 2020', 23: 'Apr 2020', 24: 'May 2020',
           25: 'Jun 2020', 26: 'Jul 2020', 27: 'Jul 2020', 28: 'Aug 2020',
           29: 'Sep 2020', 30: 'Oct 2020', 31: 'Nov 2020', 32: 'Dec 2020',
           33: 'Jan 2021', 34: 'Feb 2021', 35: 'Mar 2021', 36: 'Apr 2021',
           37: 'Apr 2021', 38: 'May 2021', 39: 'Jun 2021', 40: 'Jul 2021',
           41: 'Aug 2021', 42: 'Sep 2021', 43: 'Oct 2021', 44: 'Nov 2021',
           45: 'Dec 2021', 46: 'Dec 2021', 47: 'Jan 2022', 48: 'Feb 2022',
           49: 'Mar 2022', 50: 'Apr 2022', 51: 'May 2022', 52: 'Jun 2022',
           53: 'Jul 2022', 54: 'Aug 2022', 55: 'Sep 2022', 56: 'Sep 2022',
           57: 'Oct 2022', 58: 'Nov 2022', 59: 'Dec 2022', 60: 'Jan 2023',
           61: 'Feb 2023', 62: 'Mar 2023', 63: 'Apr 2023', 64: 'May 2023',
           65: 'Jun 2023', 66: 'Jul 2023', 67: 'Jul 2023', 68: 'Aug 2023',
           69: 'Sep 2023', 70: 'Oct 2023', 71: 'Nov 2023', 72: 'Dec 2023',
           73: 'Jan 2024', 74: 'Jan 2024', 75: 'Feb 2024', 76: 'Mar 2024',
           77: 'Apr 2024', 78: 'May 2024', 79: 'Jun 2024', 80: 'Jul 2024',
           81: 'Aug 2024', 82: 'Sep 2024', 83: 'Oct 2024'}

##################################################################

warnings.filterwarnings("ignore", category=FITSFixedWarning)

if not test and not savefig:
    warnings.warn('Not testing but not saving any output.')

if ecliptic_coords and galactic_coords:
    raise Exception("can't use both ecliptic and galactic coords.")

cnorm = colors.LogNorm(vmin=vmin, vmax=vmax)
if color == 'gray':
    # set up our custom colormap, which is a subset of the mpl map 'gray'.
    # we use truncate_colormap() to remove the blackest part of the map
    # so that even the dark areas show up against a pure black background.
    cmap = 'gray'
    # use only the latter part of the original colormap
    cmap = truncate_colormap(plt.get_cmap(cmap), minval=0.18, maxval=1.0)
elif color == 'blue':
    img = imread('bluecmap.png')
    pic = img[:, 0, :]
    cmap = colors.LinearSegmentedColormap.from_list('NASA_blue', pic[::-1, :],
                                                    N=pic.shape[0])
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
        fcam = int(os.path.split(ifile)[1].split('-')[2])
        fccd = int(os.path.split(ifile)[1].split('-')[3])
        # decide if we want to use it
        if (((fsec < 14) or ((fsec > 26) and (fsec < 40)) or (fsec > 60)) and
                hemisphere in ['both', 'south']):
            files.append(ifile)
        elif ((((fsec > 13) and (fsec < 27)) or
               ((fsec > 39) and (fsec < 42)) or
               ((fsec > 46) and (fsec < 61))) and
              hemisphere in ['both', 'north']):
            files.append(ifile)
        elif (fsec > 41) and (fsec < 47):
            if (hemisphere in ['both', 'north'] and
                    ((fcam in [1, 2] and fccd in [2, 3]) or
                     (fcam in [3, 4] and fccd in [1, 4]))):
                files.append(ifile)
            if (hemisphere in ['both', 'south'] and
                    ((fcam in [1, 2] and fccd in [1, 4]) or
                     (fcam in [3, 4] and fccd in [2, 3]))):
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
    binning = False

# remove any previous images in the gif subdirectory
if makegif:
    prev = glob(os.path.join(figdir, '*png'))
    for iprev in prev:
        os.remove(iprev)

# anything we want to test
if test:
    files = glob(os.path.join(datadir, f'*s0062-4*fits'))
    # files += glob(os.path.join(datadir, f'*s0059-2-3*fits'))
    # files += glob(os.path.join(datadir, f'*s001[56]-1-[12]*fits'))
    # files += glob(os.path.join(datadir, f'*s0019-3-[12]*fits'))

    # files += glob(os.path.join(datadir, f'*s0023-2-[12]*fits'))
    # files += glob(os.path.join(datadir, f'*s0022-1-[12]*fits'))
    # files += glob(os.path.join(datadir, f'*s0021*fits'))
    files.sort()

    kepfiles = []
    # kepfiles = glob(os.path.join(datadir, f'k*c08*fits'))
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
    fthru = int(len(lines) * frac)
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
    fscl = 95
    if hemisphere == 'both':
        xinch = 150
        yinch = 75
    else:
        xinch = 100
        yinch = 100
    fsz = int(160 * fscl / 100.)
    sfsz = int(175 * fscl / 100.)
    tfsz = int(200 * fscl / 100.)
elif fullres:
    fscl = 400
    if hemisphere == 'both':
        xinch = 600
        yinch = 300
    else:
        xinch = 400
        yinch = 400
    fsz = int(160 * fscl / 100.)
    sfsz = int(175 * fscl / 100.)
    tfsz = int(200 * fscl / 100.)
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

# the data coordinates are lat/lon in a grid
data_tr = ccrs.PlateCarree()

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
    # load the edges of all the outer CCDs and invisibly plot them
    # so that after just 1 sector, the plot isn't artificially
    # zoomed to just that one sector.
    for edgefile in edgefiles:
        elat, elon = np.loadtxt(edgefile, unpack=True)
        if not test:
            plt.scatter(elon, elat, c='w', alpha=0.01, zorder=-5,
                        marker='.', s=1, transform=data_tr)

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
    print(f'{now}. Processing Kepler/K2 image {ict + 1} of {len(kepfiles)}.')
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
            xinds = np.arange(-0.5, data.shape[0] - 0.4)
            yinds = np.arange(-0.5, data.shape[1] - 0.4)
            mesh = np.meshgrid(xinds, yinds, indexing='ij')

            # transofrm to actual sky coordinates
            lon, lat = wcs.all_pix2world(mesh[1].flatten(),
                                         mesh[0].flatten(), 0)
            lon = lon.reshape(mesh[0].shape)
            lat = lat.reshape(mesh[1].shape)

            # detect bad modules
            if (lon[0, 0] == 0.5) and (lat[0, 0] == 0.5):
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
            if galactic_coords:
                icrs = SkyCoord(ra=lon, dec=lat, frame='icrs', unit='deg')
                galactic = icrs.transform_to(Galactic)
                lon = galactic.l.value * 1
                lat = galactic.b.value * 1

            # lon must be between -180 and 180 instead of 0 to 360
            lon -= 180.
            # because in astronomy images, plots have east on the left,
            # switch east and west
            lon *= -1.

            # rebin from Kepler's 4 arcsec to 16 arcsec pixels, closer to TESS
            if binning:
                bb = 8
                data = data[:, :-4]
                lon = lon[:, :-4]
                lat = lat[:, :-4]
            else:
                bb = 4
            data = rebin(data, (data.shape[0] // bb, data.shape[1] // bb))
            lat = lat[::bb, ::bb]
            lon = lon[::bb, ::bb]

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
                lmin = (((cenlon - 178) + 180) % 360) - 180
                lmax = (((cenlon + 178) + 180) % 360) - 180
                wlon = ((lon - cenlon + 180) % 360) - 180
                if wrap and (lon.max() > lmax) and (lon.min() < lmin):
                    # find the problem areas that wrap around in longitude
                    bad = ((np.abs(wlon[:-1, :-1] - wlon[:-1, 1:]) > 355.) |
                           (np.abs(wlon[:-1, :-1] - wlon[1:, :-1]) > 355.) |
                           (np.abs(wlon[:-1, :-1] - wlon[1:, 1:]) > 355.))
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

imct = 0
# loop through every image and create the mosaic
for ii, ifile in enumerate(files):
    with fits.open(ifile) as ff:
        wcs = WCS(ff[1].header)
        data = ff[1].data * 1

        # get the coordinates of every data point in CCD coordinates
        xinds = np.arange(-0.5, data.shape[0] - 0.4)
        yinds = np.arange(-0.5, data.shape[1] - 0.4)
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
        if galactic_coords:
            icrs = SkyCoord(ra=lon, dec=lat, frame='icrs', unit='deg')
            galactic = icrs.transform_to(Galactic)
            lon = galactic.l.value * 1
            lat = galactic.b.value * 1

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
        print(
            f'{now}. Processing image {ii + 1} of {len(files)}: Sector {isec} '
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
        elif isec == 19 and icam == 1 and iccd == 2:
            xx = np.where((data[200:450, 1750:2000] >= 1.1 * vmin) &
                          (data[200:450, 1750:2000] <= 450))
            rr = np.random.rand(*data[200:450, 1750:2000].shape)
            rr *= 0.15 * vmin
            rr += vmin
            mins = rr[xx]
            data[200:450, 1750:2000][xx] = np.clip(data[200:450, 1750:2000][xx]
                                                   - 250, mins, None)
        elif isec == 19 and icam == 3 and iccd == 3:
            # a reflection in the CVZ only in S19
            data = data[360:, :]
            lat = lat[360:, :]
            lon = lon[360:, :]
        elif isec == 20 and icam == 1 and iccd == 2:
            xx = np.where(data[100:350, 1700:1950] >= vmin)
            data[100:350, 1700:1950][xx] /= 1.5
        elif isec == 21 and icam == 2 and iccd == 1:
            # a reflection in the CVZ only in S21
            data = data[226:, :]
            lat = lat[226:, :]
            lon = lon[226:, :]
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
        elif isec == 39 and icam == 3 and iccd == 2:
            data[:1000, 1300:] = np.nan
        elif isec == 39 and icam == 2 and iccd == 1:
            data[:800, :1000] = np.nan
        elif isec == 39 and icam == 1 and iccd == 2:
            data[:2000, 1850:] = np.nan
        elif isec == 39 and icam == 1 and iccd == 1:
            data[:500, :800] = np.nan
        elif isec == 39 and icam == 2 and iccd == 4:
            data[:750, 1450:] = np.nan
        elif isec == 39 and icam == 2 and iccd == 3:
            data[250:, :500] = np.nan
        elif isec == 39 and icam == 1 and iccd == 4:
            data[:750, 1200:] = np.nan
        elif isec == 40 and icam == 4 and iccd == 3:
            data[:300, :800] = np.nan
        elif isec == 40 and icam == 4 and iccd == 4:
            data[:300, 1400:] = np.nan
        elif isec == 41 and icam == 1 and iccd == 1:
            data[50:300, :100] = np.nan
        elif isec == 41 and icam == 1 and iccd == 2:
            data[:200, 1700:1900] = np.nan
        elif isec == 41 and icam == 2 and iccd == 3:
            data[:300, :300] = np.nan
        elif isec == 41 and icam == 2 and iccd == 4:
            data[175:400, 1830:] = np.nan
        elif isec == 41 and icam == 3 and iccd == 1:
            data[:300, :400] = np.nan
        elif isec == 41 and icam == 4 and iccd == 3:
            data[200:400, :100] = np.nan
        elif isec == 41 and icam == 4 and iccd == 4:
            data[:100, 1400:] = np.nan
        elif isec == 42 and icam == 3 and iccd == 3:
            data[:150, :200] = np.nan
        elif isec == 42 and icam == 4 and iccd == 4:
            data[:1500, -1000:] = np.nan
        elif isec == 43 and icam == 1 and iccd == 4:
            data[:300, 1700:] = np.nan
        elif isec == 43 and icam == 2 and iccd == 2:
            data[1800:, :400] = np.nan
        elif isec == 44 and icam == 1 and iccd == 1:
            data[:300, 200:800] = np.nan
        elif isec == 44 and icam == 2 and iccd == 1:
            data[600:1100, :400] = np.nan
        elif isec == 44 and icam == 2 and iccd == 4:
            data[:350, 1400:] = np.nan
        elif isec == 44 and icam == 3 and iccd == 1:
            data[:400, :200] = np.nan
        elif isec == 44 and icam == 3 and iccd == 2:
            data[300:1000, 1600:] = np.nan
        elif isec == 44 and icam == 3 and iccd == 3:
            data[200:400, :200] = np.nan
        elif isec == 44 and icam == 3 and iccd == 4:
            data[100:300, 1900:] = np.nan
        elif isec == 44 and icam == 4 and iccd == 2:
            data[:500, 1600:] = np.nan
        elif isec == 44 and icam == 4 and iccd == 3:
            data[:400, :700] = np.nan
        elif isec == 44 and icam == 4 and iccd == 4:
            data[:200, 1400:1800] = np.nan
        elif isec == 45 and icam == 1 and iccd == 1:
            data[:200, :300] = np.nan
        elif isec == 45 and icam == 1 and iccd == 2:
            data[:1100, 1000:] = np.nan
        elif isec == 45 and icam == 1 and iccd == 3:
            data[:400, :400] = np.nan
        elif isec == 45 and icam == 1 and iccd == 4:
            data[:800, 1700:] = np.nan
            data[:200, 1350:1550] = np.nan
        elif isec == 45 and icam == 2 and iccd == 1:
            data[:250, 150:450] = np.nan
        elif isec == 45 and icam == 2 and iccd == 2:
            data[:300, 1600:1950] = np.nan
        elif isec == 45 and icam == 2 and iccd == 3:
            data[:350, :350] = np.nan
            data[:150, 450:750] = np.nan
        elif isec == 45 and icam == 2 and iccd == 4:
            data[:200, 1350:1550] = np.nan
            data[100:800, 1800:] = np.nan
        elif isec == 45 and icam == 3 and iccd == 2:
            data[:200, 1300:1600] = np.nan
            data[200:400, 1900:] = np.nan
        elif isec == 45 and icam == 3 and iccd == 3:
            data[:200, 500:650] = np.nan
            data[:350, :150] = np.nan
        elif isec == 45 and icam == 4 and iccd == 3:
            data[:200, 500:700] = np.nan
            data[190:350, :100] = np.nan
        elif isec == 45 and icam == 4 and iccd == 4:
            data[:200, 1300:1600] = np.nan
            data[:300, 1900:] = np.nan
        elif isec == 46 and icam == 1 and iccd == 3:
            data[:200, 500:700] = np.nan
        elif isec == 46 and icam == 2 and iccd == 3:
            data[:250, :100] = np.nan
        elif isec == 46 and icam == 2 and iccd == 4:
            data[:450, 1800:] = np.nan
        elif isec == 46 and icam == 4 and iccd == 1:
            data[70:450, :200] = np.nan
        elif isec == 46 and icam == 4 and iccd == 2:
            data[:200, 1900:] = np.nan
        elif isec == 46 and icam == 4 and iccd == 3:
            data[:700, :700] = np.nan
        elif isec == 46 and icam == 4 and iccd == 4:
            data[:200, 1400:1850] = np.nan
            data[500:700, 1900:] = np.nan
        elif isec == 47 and icam == 1 and iccd == 3:
            data[:200, :600] = np.nan
        elif isec == 47 and icam == 1 and iccd == 4:
            data[:200, 1700:] = np.nan
        elif isec == 47 and icam == 3 and iccd == 2:
            data[:350, 1700:] = np.nan
        elif isec == 48 and icam == 2 and iccd == 4:
            xx = np.where(data[:160, 1570:1780] >= vmin)
            data[:160, 1570:1780][xx] /= 2.5
        elif isec == 49 and icam == 2 and iccd == 3:
            data[:200, 100:400] = np.nan
        elif isec == 50 and icam == 1 and iccd == 2:
            data[:200, 1650:1900] = np.nan
        elif isec == 50 and icam == 2 and iccd == 2:
            data[100:300, 1550:1800] = np.nan
        elif isec == 51 and icam == 3 and iccd == 3:
            data[75:200, :150] = np.nan
        elif isec == 51 and icam == 1 and iccd == 4:
            xx = np.where(data[20:315, 1850:] >= 1. * vmin)
            cclip = np.clip(data[20:315, 1850:][xx] - 100, a_min=1. * vmin,
                            a_max=np.inf)
            data[20:315, 1850:][xx] = cclip + np.random.randint(-2, 25,
                                                                xx[0].size)
            xx = np.where(data[20:50, 1970:] >= 1. * vmin)
            cclip = np.clip(data[20:50, 1970:][xx] - 200, a_min=1. * vmin,
                            a_max=np.inf)
            data[20:50, 1970:][xx] = cclip + np.random.randint(-2, 25,
                                                               xx[0].size)
        elif isec == 52 and icam == 3 and iccd == 2:
            data[50:200, 1850:] = np.nan
        elif isec == 53 and icam == 3 and iccd == 3:
            data[50:350, 30:300] = np.nan
        elif isec == 53 and icam == 4 and iccd == 4:
            data[:200, 1850:] = np.nan
        elif isec == 55 and icam == 2 and iccd == 1:
            data[:150, :200] = np.nan
        elif isec == 55 and icam == 3 and iccd == 4:
            data[:150, 1600:1900] = np.nan
        elif isec == 55 and icam == 3 and iccd == 1:
            data[:500, :700] = np.nan
        elif isec == 56 and icam == 1 and iccd == 3:
            data[:150, 50:200] = np.nan
        elif isec == 56 and icam == 1 and iccd == 4:
            data[50:150, 1900:] = np.nan
        elif isec == 56 and icam == 2 and iccd == 4:
            data[:10, 940:1120] = np.nan
            data[:50, 940:1120] -= 150.
            data[:50, 940:1120] = np.clip(data[:50, 940:1120], a_min=50.,
                                          a_max=None)
        elif isec == 58 and icam == 1 and iccd == 4:
            data[400:1100, 800:] = np.nan
            data[:200, 1800:] = np.nan
        elif isec == 58 and icam == 2 and iccd == 2:
            data[:500, 1400:] = np.nan
        elif isec == 58 and icam == 2 and iccd == 4:
            data[:250, 1400:] = np.nan
        elif isec == 58 and icam == 3 and iccd == 1:
            data[:400, :800] = np.nan
        elif isec == 59 and icam == 1 and iccd == 3:
            data[:600, :600] = np.nan
        elif isec == 59 and icam == 3 and iccd == 3:
            data[250:600, :200] = np.nan
        elif isec == 61 and icam == 1 and iccd == 3:
            data[:200, :300] = np.nan
        elif isec == 61 and icam == 2 and iccd == 2:
            data[:400, 1500:] = np.nan
        elif isec == 61 and icam == 2 and iccd == 4:
            data[:600, 1600:] = np.nan
        elif isec == 61 and icam == 3 and iccd == 2:
            data[:400, 1550:1950] = np.nan

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
        if isec == 39 and icam == 3 and iccd in [2, 3]:
            data += 10
        # galactic center is too bright
        if isec in [12, 39] and icam == 1:
            data -= 20
        if isec == 13 and icam == 1 and iccd in [1, 4]:
            data -= 30
        if isec == 40 and icam == 1 and iccd == 3:
            data[:200, 250:550] -= 40
        if isec == 41 and icam == 1 and iccd == 4:
            data[:200, 1600:1950] -= 40
        if isec == 42 and icam == 1 and iccd == 2:
            data[:200, 1900:] -= 80
        if isec == 42 and icam == 4 and iccd in [2, 3]:
            data -= 25
        if isec == 42 and icam == 4 and iccd in [1, 4]:
            data += 15
        if isec == 42 and icam == 1 and iccd in [1, 2, 3]:
            data -= 30
        if isec == 53 and icam == 4 and iccd in [3, 4]:
            data -= 20
        if isec == 42 and icam == 1 and iccd in [4]:
            # top, bottom
            data += np.ones_like(data) * np.linspace(-15, 5, data.shape[0])
            # right, left
            data += (np.ones_like(data) * np.linspace(20, -5, data.shape[0])).T

        if isec == 43 and icam == 3 and iccd in [1]:
            # right, left
            data += (np.ones_like(data) * np.linspace(0, -15, data.shape[0])).T
            # top, bottom
            data += np.ones_like(data) * np.linspace(40, 0, data.shape[0])

        if isec == 43 and icam == 3 and iccd in [4]:
            # left, right
            data += (np.ones_like(data) * np.linspace(-10, -5, data.shape[0])).T
            # bottom, top
            data += np.ones_like(data) * np.linspace(-15, 20, data.shape[0])
        if isec == 43 and icam == 4 and iccd in [1]:
            # right, left
            data += (np.ones_like(data) * np.linspace(-5, -50, data.shape[0])).T
            # top, bottom
            data += np.ones_like(data) * np.linspace(20, -15, data.shape[0])

        if isec == 43 and icam == 4 and iccd in [4]:
            data -= 30
            data += np.ones_like(data) * np.linspace(-30, 0, data.shape[0])
            data += (np.ones_like(data) * np.linspace(-15, 0, data.shape[0])).T
        if isec == 43 and icam == 3 and iccd in [2]:
            # right, left
            data += (np.ones_like(data) * np.linspace(0, -8, data.shape[0])).T
            # top, bottom
            data += np.ones_like(data) * np.linspace(0, 5, data.shape[0])
        if isec == 43 and icam == 3 and iccd in [3]:
            # bottom, top
            data += np.ones_like(data) * np.linspace(7, -20, data.shape[0])
            # left, right
            data += (np.ones_like(data) * np.linspace(-15, -5, data.shape[0])).T
        if isec == 43 and icam == 4 and iccd in [2]:
            data -= 20
            # top, bottom
            data += np.ones_like(data) * np.linspace(-10, 5, data.shape[0])
            # right, left
            data += (np.ones_like(data) * np.linspace(0, -25, data.shape[0])).T
        if isec == 43 and icam == 4 and iccd in [3]:
            data -= 40
            # bottom, top
            data += np.ones_like(data) * np.linspace(5, -20, data.shape[0])
            # left, right
            data += (np.ones_like(data) * np.linspace(-40, -5, data.shape[0])).T

        if isec == 44 and icam == 1 and iccd == 1:
            # left, right
            data += (np.ones_like(data) * np.linspace(5, -25, data.shape[0])).T
            # bottom, top
            data += np.ones_like(data) * np.linspace(15, -10, data.shape[0])
        if isec == 44 and icam == 1 and iccd == 2:
            # left, right
            data += (np.ones_like(data) * np.linspace(0, -15, data.shape[0])).T
            # bottom, top
            data += np.ones_like(data) * np.linspace(-10, 30, data.shape[0])
        if isec == 44 and icam == 1 and iccd == 3:
            # top, bottom
            data += np.ones_like(data) * np.linspace(20, -30, data.shape[0])
            # right, left
            data += (np.ones_like(data) * np.linspace(0, -5, data.shape[0])).T
        if isec == 44 and icam == 1 and iccd == 4:
            # right, left
            data += (np.ones_like(data) * np.linspace(-5, -10, data.shape[0])).T
            # top, bottom
            data += np.ones_like(data) * np.linspace(-25, 0, data.shape[0])

        if isec == 44 and icam == 2 and iccd == 1:
            # left, right
            data += (np.ones_like(data) * np.linspace(-15, -10,
                                                      data.shape[0])).T
            # bottom, top
            data += np.ones_like(data) * np.linspace(19, 0, data.shape[0])
        if isec == 44 and icam == 2 and iccd == 2:
            # left, right
            data += (np.ones_like(data) * np.linspace(-5, 0, data.shape[0])).T
            # bottom, top
            data += np.ones_like(data) * np.linspace(-10, 25, data.shape[0])
        if isec == 44 and icam == 2 and iccd == 3:
            # top, bottom
            data += np.ones_like(data) * np.linspace(30, -10, data.shape[0])
            # right, left
            data += (np.ones_like(data) * np.linspace(-5, 0, data.shape[0])).T
        if isec == 44 and icam == 2 and iccd == 4:
            # right, left
            data += (np.ones_like(data) * np.linspace(0, -10, data.shape[0])).T
            # top, bottom
            data += np.ones_like(data) * np.linspace(0, 26, data.shape[0])

        if isec == 44 and icam == 3 and iccd == 1:
            # right, left
            data += (np.ones_like(data) * np.linspace(0, -15, data.shape[0])).T
            # top, bottom
            data += np.ones_like(data) * np.linspace(30, 0, data.shape[0])
        if isec == 44 and icam == 3 and iccd == 2:
            # top, bottom
            data += np.ones_like(data) * np.linspace(-15, 0, data.shape[0])
            # right, left
            data += (np.ones_like(data) * np.linspace(5, 5, data.shape[0])).T
        if isec == 44 and icam == 3 and iccd == 3:
            # bottom, top
            data += np.ones_like(data) * np.linspace(0, -20, data.shape[0])
            # left, right
            data += (np.ones_like(data) * np.linspace(-5, 0, data.shape[0])).T
        if isec == 44 and icam == 3 and iccd == 4:
            # left, right
            data += (np.ones_like(data) * np.linspace(-10, -5, data.shape[0])).T
            # bottom, top
            data += np.ones_like(data) * np.linspace(-10, 20, data.shape[0])

        if isec == 44 and icam == 4 and iccd == 1:
            # right, left
            data += (np.ones_like(data) * np.linspace(0, -20, data.shape[0])).T
            # top, bottom
            data += np.ones_like(data) * np.linspace(0, -30, data.shape[0])
        if isec == 44 and icam == 4 and iccd == 2:
            # top, bottom
            data += np.ones_like(data) * np.linspace(-30, 10, data.shape[0])
            # right, left
            data += (np.ones_like(data) * np.linspace(0, -20, data.shape[0])).T
        if isec == 44 and icam == 4 and iccd == 3:
            # bottom, top
            data += np.ones_like(data) * np.linspace(-10, -50, data.shape[0])
            # left, right
            data += (np.ones_like(data) * np.linspace(-10, 0, data.shape[0])).T
        if isec == 44 and icam == 4 and iccd == 4:
            # left, right
            data += (np.ones_like(data) * np.linspace(-20, 0, data.shape[0])).T
            # bottom, top
            data += np.ones_like(data) * np.linspace(-60, 0, data.shape[0])

        if isec == 45 and icam == 1 and iccd == 1:
            # left, right
            data += (np.ones_like(data) * np.linspace(10, -10, data.shape[0])).T
            # bottom, top
            data += np.ones_like(data) * np.linspace(20, -10, data.shape[0])
        if isec == 45 and icam == 1 and iccd == 2:
            # left, right
            data += (np.ones_like(data) * np.linspace(10, -10, data.shape[0])).T
            # bottom, top
            data += np.ones_like(data) * np.linspace(0, 20, data.shape[0])
        if isec == 45 and icam == 1 and iccd == 3:
            # top, bottom
            data += np.ones_like(data) * np.linspace(10, 0, data.shape[0])
            # right, left
            data += (np.ones_like(data) * np.linspace(0, -10, data.shape[0])).T
        if isec == 45 and icam == 1 and iccd == 4:
            # right, left
            data += (np.ones_like(data) * np.linspace(10, -10, data.shape[0])).T
            # top, bottom
            data += np.ones_like(data) * np.linspace(-20, 10, data.shape[0])

        if isec == 45 and icam == 2 and iccd == 1:
            # left, right
            data += (np.ones_like(data) * np.linspace(0, -10, data.shape[0])).T
            # bottom, top
            data += np.ones_like(data) * np.linspace(0, -10, data.shape[0])
        if isec == 45 and icam == 2 and iccd == 2:
            # left, right
            data += (np.ones_like(data) * np.linspace(0, -10, data.shape[0])).T
            # bottom, top
            data += np.ones_like(data) * np.linspace(-10, 0, data.shape[0])
        if isec == 45 and icam == 2 and iccd == 3:
            # top, bottom
            data += np.ones_like(data) * np.linspace(20, 0, data.shape[0])
            # right, left
            data += (np.ones_like(data) * np.linspace(10, -20, data.shape[0])).T
        if isec == 45 and icam == 2 and iccd == 4:
            # right, left
            data += (np.ones_like(data) * np.linspace(10, -10, data.shape[0])).T
            # top, bottom
            data += np.ones_like(data) * np.linspace(0, 20, data.shape[0])

        if isec == 45 and icam == 3 and iccd == 1:
            # right, left
            data += (np.ones_like(data) * np.linspace(0, -10, data.shape[0])).T
            # top, bottom
            data += np.ones_like(data) * np.linspace(0, -10, data.shape[0])
        if isec == 45 and icam == 3 and iccd == 2:
            # top, bottom
            data += np.ones_like(data) * np.linspace(-10, 0, data.shape[0])
            # right, left
            data += (np.ones_like(data) * np.linspace(0, -10, data.shape[0])).T
        if isec == 45 and icam == 3 and iccd == 3:
            # bottom, top
            data += np.ones_like(data) * np.linspace(0, -40, data.shape[0])
            # left, right
            data += (np.ones_like(data) * np.linspace(0, 0, data.shape[0])).T
        if isec == 45 and icam == 3 and iccd == 4:
            # left, right
            data += (np.ones_like(data) * np.linspace(0, 0, data.shape[0])).T
            # bottom, top
            data += np.ones_like(data) * np.linspace(-40, 0, data.shape[0])

        if isec == 45 and icam == 4 and iccd == 1:
            # right, left
            data += (np.ones_like(data) * np.linspace(0, -50, data.shape[0])).T
            # top, bottom
            data += np.ones_like(data) * np.linspace(0, -40, data.shape[0])
        if isec == 45 and icam == 4 and iccd == 2:
            # top, bottom
            data += np.ones_like(data) * np.linspace(-40, 0, data.shape[0])
            # right, left
            data += (np.ones_like(data) * np.linspace(0, -50, data.shape[0])).T
        if isec == 45 and icam == 4 and iccd == 3:
            # bottom, top
            data += np.ones_like(data) * np.linspace(-20, -90, data.shape[0])
            # left, right
            data += (np.ones_like(data) *
                     np.linspace(-40, -20, data.shape[0])).T
        if isec == 45 and icam == 4 and iccd == 4:
            data -= 50
            # left, right
            data += (np.ones_like(data) *
                     np.linspace(-20, -10, data.shape[0])).T
            # bottom, top
            data += np.ones_like(data) * np.linspace(-60, 0, data.shape[0])

        if isec == 46 and icam == 1 and iccd == 1:
            # left, right
            data += (np.ones_like(data) * np.linspace(0, -10, data.shape[0])).T
            # bottom, top
            data += np.ones_like(data) * np.linspace(10, -10, data.shape[0])
        if isec == 46 and icam == 1 and iccd == 2:
            # left, right
            data += (np.ones_like(data) * np.linspace(0, -10, data.shape[0])).T
            # bottom, top
            data += np.ones_like(data) * np.linspace(-10, 0, data.shape[0])
        if isec == 46 and icam == 1 and iccd == 3:
            # top, bottom
            data += np.ones_like(data) * np.linspace(0, -10, data.shape[0])
            # right, left
            data += (np.ones_like(data) * np.linspace(0, -10, data.shape[0])).T
        if isec == 46 and icam == 1 and iccd == 4:
            # right, left
            data += (np.ones_like(data) * np.linspace(10, -10, data.shape[0])).T
            # top, bottom
            data += np.ones_like(data) * np.linspace(-10, 0, data.shape[0])

        if isec == 46 and icam == 3 and iccd == 1:
            # right, left
            data += (np.ones_like(data) * np.linspace(0, -20, data.shape[0])).T
            # top, bottom
            data += np.ones_like(data) * np.linspace(0, -10, data.shape[0])
        if isec == 46 and icam == 3 and iccd == 2:
            # top, bottom
            data += np.ones_like(data) * np.linspace(-10, 0, data.shape[0])
            # right, left
            data += (np.ones_like(data) * np.linspace(0, -10, data.shape[0])).T
        if isec == 46 and icam == 3 and iccd == 3:
            # bottom, top
            data += np.ones_like(data) * np.linspace(0, -30, data.shape[0])
            # left, right
            data += (np.ones_like(data) * np.linspace(0, 0, data.shape[0])).T
        if isec == 46 and icam == 3 and iccd == 4:
            # left, right
            data += (np.ones_like(data) * np.linspace(0, 0, data.shape[0])).T
            # bottom, top
            data += np.ones_like(data) * np.linspace(-30, 0, data.shape[0])

        if isec == 46 and icam == 4 and iccd == 1:
            # right, left
            data += (np.ones_like(data) * np.linspace(0, -40, data.shape[0])).T
            # top, bottom
            data += np.ones_like(data) * np.linspace(0, -40, data.shape[0])
        if isec == 46 and icam == 4 and iccd == 2:
            # top, bottom
            data += np.ones_like(data) * np.linspace(-40, 0, data.shape[0])
            # right, left
            data += (np.ones_like(data) * np.linspace(0, -40, data.shape[0])).T
        if isec == 46 and icam == 4 and iccd == 3:
            data -= 20
            # bottom, top
            data += np.ones_like(data) * np.linspace(0, -70, data.shape[0])
            # left, right
            data += (np.ones_like(data) * np.linspace(-30, 0, data.shape[0])).T
        if isec == 46 and icam == 4 and iccd == 4:
            data -= 40
            # left, right
            data += (np.ones_like(data) * np.linspace(-30, 0, data.shape[0])).T
            # bottom, top
            data += np.ones_like(data) * np.linspace(-60, 0, data.shape[0])

        if isec == 47 and icam == 1 and iccd == 1:
            # top, bottom
            data += (np.ones_like(data) * np.linspace(0, 0, data.shape[0])).T
            # right, left
            data += np.ones_like(data) * np.linspace(0, 0, data.shape[0])
        if isec == 47 and icam == 1 and iccd == 2:
            # right, left
            data += (np.ones_like(data) * np.linspace(10, 0, data.shape[0])).T
            # bottom, top
            data += np.ones_like(data) * np.linspace(0, 0, data.shape[0])
        if isec == 47 and icam == 1 and iccd == 3:
            # right, left
            data += np.ones_like(data) * np.linspace(20, 0, data.shape[0])
            # bottom, top
            data += (np.ones_like(data) * np.linspace(0, 0, data.shape[0])).T
        if isec == 47 and icam == 1 and iccd == 4:
            # bottom, top
            data += (np.ones_like(data) * np.linspace(0, 0, data.shape[0])).T
            # right, left
            data += np.ones_like(data) * np.linspace(0, 0, data.shape[0])

        if isec == 52 and icam == 4 and iccd == 1:
            # bottom, top
            data += (np.ones_like(data) * np.linspace(0, 0, data.shape[0])).T
            # right, left
            data += np.ones_like(data) * np.linspace(0, 0, data.shape[0])
        if isec == 52 and icam == 4 and iccd == 2:
            # bottom, top
            data += (np.ones_like(data) * np.linspace(0, -20, data.shape[0])).T
            # right, left
            data += np.ones_like(data) * np.linspace(0, -10, data.shape[0])
        if isec == 52 and icam == 4 and iccd == 3:
            # left, right
            data += np.ones_like(data) * np.linspace(-30, -10, data.shape[0])
            # top, bottom
            data += (np.ones_like(data) * np.linspace(-50, 0, data.shape[0])).T
        if isec == 52 and icam == 4 and iccd == 4:
            # top, bottom
            data += (np.ones_like(data) * np.linspace(-40, 0, data.shape[0])).T
            # left, right
            data += np.ones_like(data) * np.linspace(-30, 0, data.shape[0])

        if corner_glow_plot:
            plt.figure()
            plt.imshow(data, norm=cnorm, cmap=cmap)

        if makefig:
            if binning:
                data = rebin(data, (data.shape[0] // 2, data.shape[1] // 2))
                lon = lon[::2, ::2]
                lat = lat[::2, ::2]
            if ii == 0:
                sectxt = f'Sector {isec}\n{secstarts[isec]}' \
                         f'\u2013{secends[isec]}'
                if printdate:
                    text = plt.text(0.98, 0.02, sectxt,
                                    transform=fig.transFigure, ha='right',
                                    va='bottom', multialignment='right',
                                    fontsize=sfsz, fontname='Carlito')
                ssec = isec
            # for wraparounds:
            lmin = (((cenlon - 178) + 180) % 360) - 180
            lmax = (((cenlon + 178) + 180) % 360) - 180
            wlon = ((lon - cenlon + 180) % 360) - 180
            if wrap and (lon.max() > lmax) and (lon.min() < lmin):
                # find the problem areas that wrap around in longitude
                bad = ((np.abs(wlon[:-1, :-1] - wlon[:-1, 1:]) > 355.) |
                       (np.abs(wlon[:-1, :-1] - wlon[1:, :-1]) > 355.) |
                       (np.abs(wlon[:-1, :-1] - wlon[1:, 1:]) > 355.))
                # mask them and just don't plot these pixels
                maskeddata = np.ma.masked_where(bad, data)
                plt.pcolormesh(lon, lat, maskeddata, norm=cnorm, alpha=1,
                               transform=data_tr, cmap=cmap)

            else:
                # plot the actual image from this CCD
                plt.pcolormesh(lon, lat, data, norm=cnorm, alpha=1,
                               transform=data_tr, cmap=cmap)

        psec = 0
        if ii > 0:
            psec = int(files[ii-1].split('-s0')[1][:3])
        if isec != psec and ii > 0:
            newsec = True
        else:
            newsec = False
        # if we're starting to plot a new sector, update the date
        if newsec and printdate and makefig:
            text.remove()
            if hemisphere == 'both' or isec < 27:
                sectxt = f'Sectors {ssec}\u2013{isec}\n{secstarts[ssec]}' \
                         f'\u2013{secends[isec]}'
            elif hemisphere == 'south':
                if isec == 27:
                    sectxt = f'Sectors 1\u201313; 27\n{secstarts[1]}' \
                             f'\u2013{secends[13]}\n{secstarts[27]}' \
                             f'\u2013{secends[27]}'
                elif isec < 40:
                    sectxt = f'Sectors 1\u201313; 27\u2013{isec}\n' \
                             f'{secstarts[1]}\u2013{secends[13]}\n' \
                             f'{secstarts[27]}\u2013{secends[isec]}'
                elif isec == 42:
                    sectxt = f'Sectors 1\u201313; 27\u201339; 42\n' \
                             f'{secstarts[1]}\u2013{secends[13]}\n' \
                             f'{secstarts[27]}\u2013{secends[27]}\n' \
                             f'{secstarts[42]}\u2013{secends[42]}'
                elif isec <= 46:
                    sectxt = f'Sectors 1\u201313; 27\u201339; 42\u2013{isec}\n'\
                             f'{secstarts[1]}\u2013{secends[13]}\n' \
                             f'{secstarts[27]}\u2013{secends[39]}\n' \
                             f'{secstarts[42]}\u2013{secends[isec]}'
                elif isec == 61:
                    sectxt = f'Sectors 1\u201313; 27\u201339;\n' \
                             f'42\u201346; {isec}\n' \
                             f'{secstarts[1]}\u2013{secends[13]}\n' \
                             f'{secstarts[27]}\u2013{secends[39]}\n' \
                             f'{secstarts[42]}\u2013{secends[46]}\n' \
                             f'{secstarts[61]}\u2013{secends[61]}'
                else:
                    sectxt = f'Sectors 1\u201313; 27\u201339;\n' \
                             f'42\u201346; 61\u2013{isec}\n' \
                             f'{secstarts[1]}\u2013{secends[13]}\n' \
                             f'{secstarts[27]}\u2013{secends[39]}\n' \
                             f'{secstarts[42]}\u2013{secends[46]}\n' \
                             f'{secstarts[61]}\u2013{secends[isec]}'
            elif hemisphere == 'north':
                if isec == 40:
                    sectxt = f'Sectors 14\u201326; 40\n{secstarts[14]}' \
                             f'\u2013{secends[26]}\n{secstarts[40]}' \
                             f'\u2013{secends[40]}'
                else:
                    sectxt = f'Sectors 14\u201326; 40\u2013{isec}\n' \
                             f'{secstarts[14]}\u2013{secends[26]}\n' \
                             f'{secstarts[40]}\u2013{secends[isec]}'
            text = plt.text(0.98, 0.02, sectxt, transform=fig.transFigure,
                            ha='right', va='bottom', multialignment='right',
                            fontsize=sfsz, fontname='Carlito')
        nsec = 0
        if ii < len(files)-1:
            nsec = int(files[ii+1].split('-s0')[1][:3])
        if isec != nsec or ii == len(files) - 1:
            endsec = True
        else:
            endsec = False
        # save the plot after each sector for the gif
        if makegif and savefig and makefig and endsec:
            if transparent:
                outfig = os.path.join(figdir,
                                      f'transp_img{imct:04d}.png')
            else:
                outfig = os.path.join(figdir, f'img{imct:04d}.png')
            plt.savefig(outfig, transparent=transparent)
            imct += 1

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
    fig = plt.figure()
    ax = plt.axes(projection='3d')

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
        imin = np.median(dat[3 * dx // 4:, 3 * dy // 4:])
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
