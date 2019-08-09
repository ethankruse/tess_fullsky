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
from scipy.stats import median_absolute_deviation as mad

datadir = os.path.join(os.path.split(__file__)[0], 'data')
figdir = os.path.join(os.path.split(__file__)[0], 'figs')

#files = glob(os.path.join(datadir, '*s0013*fits'))
#files = glob(os.path.join(datadir, '*s0011-3*fits'))
#files = glob(os.path.join(datadir, '*4-4*fits'))
#files = glob(os.path.join(datadir, '*-1-3*fits')) + glob(os.path.join(datadir, '*-1-4*fits'))
#files = glob(os.path.join(datadir, '*-1-4*fits'))
files = glob(os.path.join(datadir, '*fits'))
files.sort()


# central longitude of the projection
cenlon = 0.
# if the projection has a dividing line where we have to split
mustsplit = False #89.5 W, 66.2 S
tr = ccrs.Orthographic(central_longitude=-89.5, central_latitude=-66.2)
#tr = ccrs.Mollweide()

# minimum and maximum flux for the colorbar
vmin = 150
vmax = 1001.
cnorm = colors.LogNorm(vmin=vmin, vmax=vmax)
cmap = 'gray'

doclean = True
cleanplot = False
badcornerfile = os.path.join(os.path.split(__file__)[0], 'bad_corners.txt')
noisefile = os.path.join(os.path.split(__file__)[0], 'noise.txt')
edgefile = os.path.join(os.path.split(__file__)[0], 'edges.txt')

test = False
makefig = True
highres = True
savefig = True
makegif = True
transparent = True
if transparent:
    fname = 'transp_ortho.png'
else:
    fname = 'ortho.png'
savefile = os.path.join(figdir, fname)

credit = 'By Ethan Kruse\n@ethan_kruse'
title = "NASA TESS's View\nof the Southern\nHemisphere"

secstarts = {1: 'Jul 2018', 2: 'Aug 2018', 3: 'Sep 2018', 4: 'Oct 2018', 
             5: 'Nov 2018', 6: 'Dec 2018', 7: 'Jan 2019', 8: 'Feb 2019',
             9: 'Feb 2019', 10: 'Mar 2019', 11: 'Apr 2019', 12: 'May 2019',
             13: 'Jun 2019'}
secends = {1: 'Aug 2018', 2: 'Sep 2018', 3: 'Oct 2018', 4: 'Nov 2018', 
           5: 'Dec 2018', 6: 'Jan 2019', 7: 'Feb 2019', 8: 'Feb 2019',
           9: 'Mar 2019', 10: 'Apr 2019', 11: 'May 2019', 12: 'Jun 2019',
           13: 'Jul 2019'}
##################################################################

if makegif:
    figdir = os.path.join(figdir, 'gif')
    prev = glob(os.path.join(figdir, '*png'))
    for iprev in prev:
        os.remove(iprev)
    
if not os.path.exists(figdir):
    os.makedirs(figdir, exist_ok=True)

if test:
    # XXX: for testing
    files = [files[0]]
    # files = files[:1]
    pass

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


def clean(data, corners=[], cleanplot=False, info=None):
    import scipy.ndimage
    from astropy.modeling import models, fitting
    import warnings
    
    # a smoother copy of the data
    fdata = data * 1
    # remove spiky stars and bad data
    fdata[fdata > 1000] = np.median(fdata)
    fdata[fdata < 50] = 50
    # smooth things out
    fdata = scipy.ndimage.uniform_filter(fdata, size=201)
    
    xinds = np.arange(0.5, data.shape[0]-0.4)
    yinds = np.arange(0.5, data.shape[1]-0.4)
    lon, lat = np.meshgrid(xinds, yinds, indexing='ij')
    
    if info is not None:
        sec, cam, ccd = info
    else:
        sec, cam, ccd = 0, 0, 0
    # make a diagnostic plot to judge bad corners
    if cleanplot and len(corners) == 0:
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot_surface(lat, lon, fdata, cmap='gray',
                           linewidth=0, antialiased=False, alpha=0.2)
        ax.text(lat[0,0],lon[0,0],50,'C1',color='red')
        ax.text(lat[-1,0],lon[-1,0],50,'C2',color='red')
        ax.text(lat[0,-1],lon[0,-1],50,'C3',color='red')
        ax.text(lat[-1,-1],lon[-1,-1],50,'C4',color='red')
        plt.title(f'Sec {sec}, Cam {cam}, CCD {ccd}')
    # fix the desired corners
    for corner in corners:
        p_init = models.Polynomial2D(degree=3)
        fit_p = fitting.LevMarLSQFitter()
        
        bp_init = models.Polynomial2D(degree=1)
        bfit_p = fitting.LevMarLSQFitter()
        
        xl, yl = fdata.shape
        
        # pick the right corner and get views of the data
        if corner==1:
            xs = lat[:xl//4,:yl//4]
            ys = lon[:xl//4,:yl//4]
            dat = fdata[:xl//4,:yl//4]
            rdat = data[:xl//4,:yl//4]
            
            bxs1 = lat[xl//4:xl*3//8,:yl*3//8]
            bys1 = lon[xl//4:xl*3//8,:yl*3//8]
            bdat1 = fdata[xl//4:xl*3//8,:yl*3//8]
            bxs2 = lat[:xl//4, yl//4:yl*3//8]
            bys2 = lon[:xl//4, yl//4:yl*3//8]
            bdat2 = fdata[:xl//4, yl//4:yl*3//8]
            bxs = np.concatenate((bxs1.flatten(), bxs2.flatten()))
            bys = np.concatenate((bys1.flatten(), bys2.flatten()))
            bdat = np.concatenate((bdat1.flatten(), bdat2.flatten()))
            
        elif corner==2:
            xs = lat[-xl//4:,:yl//4]
            ys = lon[-xl//4:,:yl//4]
            dat = fdata[-xl//4:,:yl//4]
            rdat = data[-xl//4:,:yl//4]
            
            bxs1 = lat[-xl*3//8:-xl//4,:yl*3//8]
            bys1 = lon[-xl*3//8:-xl//4,:yl*3//8]
            bdat1 = fdata[-xl*3//8:-xl//4,:yl*3//8]
            bxs2 = lat[-xl//4:, yl//4:yl*3//8]
            bys2 = lon[-xl//4:, yl//4:yl*3//8]
            bdat2 = fdata[-xl//4:, yl//4:yl*3//8]
            bxs = np.concatenate((bxs1.flatten(), bxs2.flatten()))
            bys = np.concatenate((bys1.flatten(), bys2.flatten()))
            bdat = np.concatenate((bdat1.flatten(), bdat2.flatten()))
        elif corner==3:
            xs = lat[:xl//4,-yl//4:]
            ys = lon[:xl//4,-yl//4:]
            dat = fdata[:xl//4,-yl//4:]
            rdat = data[:xl//4,-yl//4:]
            
            bxs1 = lat[xl//4:xl*3//8,-yl*3//8:]
            bys1 = lon[xl//4:xl*3//8,-yl*3//8:]
            bdat1 = fdata[xl//4:xl*3//8,-yl*3//8:]
            bxs2 = lat[:xl//4, -yl*3//8:-yl//4]
            bys2 = lon[:xl//4, -yl*3//8:-yl//4]
            bdat2 = fdata[:xl//4, -yl*3//8:-yl//4]
            bxs = np.concatenate((bxs1.flatten(), bxs2.flatten()))
            bys = np.concatenate((bys1.flatten(), bys2.flatten()))
            bdat = np.concatenate((bdat1.flatten(), bdat2.flatten()))
        elif corner==4:
            xs = lat[-xl//4:,-yl//4:]
            ys = lon[-xl//4:,-yl//4:]
            dat = fdata[-xl//4:,-yl//4:]
            rdat = data[-xl//4:,-yl//4:]
            
            bxs1 = lat[-xl*3//8:-xl//4,-yl*3//8:]
            bys1 = lon[-xl*3//8:-xl//4,-yl*3//8:]
            bdat1 = fdata[-xl*3//8:-xl//4,-yl*3//8:]
            bxs2 = lat[-xl//4:, -yl*3//8:-yl//4]
            bys2 = lon[-xl//4:, -yl*3//8:-yl//4]
            bdat2 = fdata[-xl//4:, -yl*3//8:-yl//4]
            bxs = np.concatenate((bxs1.flatten(), bxs2.flatten()))
            bys = np.concatenate((bys1.flatten(), bys2.flatten()))
            bdat = np.concatenate((bdat1.flatten(), bdat2.flatten()))
        else:
            raise Exception('bad corner')

        # Ignore model linearity warning from the fitter
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            # run the fit
            poly = fit_p(p_init, xs, ys, dat)
            bpoly = bfit_p(bp_init, bxs, bys, bdat)

        mod = poly(xs, ys)
        bg = bpoly(xs, ys)
        
        # diagnostic plots to make sure it's working
        if cleanplot:
            from mpl_toolkits.mplot3d import Axes3D
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            ax.plot_surface(lat, lon, fdata, cmap='gray',
                               linewidth=0, antialiased=False, alpha=0.3)
            
            ax.plot_surface(xs, ys, mod, cmap='viridis',
                            linewidth=0, antialiased=False, alpha=0.2)
        
            dat -= mod
            dat += bg
            
            ax.plot_surface(lat, lon, fdata, cmap='gray',
                               linewidth=0, antialiased=False, alpha=0.6)
            
            ax.text(lat[0,0],lon[0,0],50,'C1',color='red')
            ax.text(lat[-1,0],lon[-1,0],50,'C2',color='red')
            ax.text(lat[0,-1],lon[0,-1],50,'C3',color='red')
            ax.text(lat[-1,-1],lon[-1,-1],50,'C4',color='red')
            plt.title(f'Sec {sec}, Cam {cam}, CCD {ccd}, Corner {corner}')
        
        # remove the trend from the actual data
        rdat -= mod
        rdat += bg
        
    return data


plt.close('all')

if highres:
    fsz = 160
    sfsz = 175
    tfsz = 200
else:
    fsz = 12
    sfsz = 13
    tfsz = 15
bsec, bcam, bccd, bcor = np.loadtxt(badcornerfile, unpack=True, ndmin=2, 
                                    delimiter=',', dtype=int)

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
        
        icam = ff[1].header['camera']
        iccd = ff[1].header['ccd']
        isec = int(ifile.split('-s0')[1][:3])
        
        imed = np.median(data)
        imad = mad(data, axis=None)
        
        print(f'{ii+1} of {len(files)}: {lon.min():.2f}, {lon.max():.2f}, {lat.min():.2f}, {lat.max():.2f}')
        print(f'median {imed:.2f}, mad {imad:.2f}, 3sigma {imed+3*imad:.2f}, 5sigma {imed+5*imad:.2f}')
        
        noises = np.loadtxt(noisefile, unpack=True, ndmin=2, delimiter=',')
        exists = np.where((noises[0,:]==isec) & (noises[1,:]==icam) & (noises[2,:]==iccd))[0]
        if len(exists) == 0:
            noises = np.concatenate((noises.T, [[isec, icam, iccd, imed, imad, imed+3.*imad, imed+5.*imad]]))
            np.savetxt(noisefile, noises, delimiter=', ', fmt='%.2f')
            
        if doclean:
            tocor = np.where((bsec==isec) & (bcam==icam) & (bccd==iccd))[0]
            print(f'Correcting corners {bcor[tocor]}')
            data = clean(data, cleanplot=cleanplot, corners=bcor[tocor],
                         info=(isec, icam, iccd))
        
        if makefig:
            if ii == 0 or test:
                if highres:
                    fig = plt.figure(figsize=(100,100))
                else:
                    fig = plt.figure()
                ax = plt.axes([0.01, 0.01, 0.99, 0.99], projection=tr)
                data_tr = ccrs.PlateCarree()

                elat, elon = np.loadtxt(edgefile, unpack=True)
                plt.scatter(elon, elat, c='w', alpha=0.01, zorder=-5, marker='.', s=1, transform=data_tr)
                
                plt.text(0.02, 0.02, credit, transform=fig.transFigure, ha='left', 
                     va='bottom', multialignment='left', fontsize=fsz, fontname='Carlito')

                plt.text(0.02, 0.98, title, transform=fig.transFigure, ha='left', 
                     va='top', multialignment='left', fontsize=tfsz, fontname='Carlito')
                sectxt = f'Sector {isec}\n{secstarts[isec]}-{secends[isec]}'
                text = plt.text(0.98, 0.02, sectxt, transform=fig.transFigure, ha='right', 
                     va='bottom', multialignment='right', fontsize=sfsz, fontname='Carlito')
                ssec = isec
            # for wraparounds:
            if lon.max() > cenlon + 120 and lon.min() < cenlon - 120 and mustsplit:
                left = np.where(lon > cenlon + 120)
                lonleft = lon * 1
                lonleft[left] = cenlon - 180.
                plt.pcolormesh(lonleft, lat, data, norm=cnorm, alpha=1, transform=data_tr, cmap=cmap)
                
                right = np.where(lon < cenlon - 120)
                lonright = lon * 1
                lonright[right] = cenlon + 180.
                plt.pcolormesh(lonright, lat, data, norm=cnorm, alpha=1, transform=data_tr, cmap=cmap)
            else:
                plt.pcolormesh(lon, lat, data, norm=cnorm, alpha=1, transform=data_tr, cmap=cmap)
            #plt.text(np.median(lon), np.median(lat), '{0}'.format(ii), transform=data_tr)
            
            
            #if test:
            #    plt.colorbar()

        if ((ii)%16) == 0 and ii > 0:
            text.remove()
            sectxt = f'Sectors {ssec}-{isec}\n{secstarts[ssec]}-{secends[isec]}'
            text = plt.text(0.98, 0.02, sectxt, transform=fig.transFigure, ha='right', 
                     va='bottom', multialignment='right', fontsize=sfsz, fontname='Carlito')
        if makegif and savefig and ii > 0 and ((ii+1)%16) == 0:
            #elat, elon = np.loadtxt(edgefile, unpack=True)
            #plt.scatter(elon, elat, c='w', alpha=0.01, zorder=-5, marker='.', s=1, transform=data_tr)
            if transparent:
                outfig = os.path.join(figdir, f'transp_img{(ii+1)//16:04d}_.png')
            else:
                outfig = os.path.join(figdir, f'img{(ii+1)//16:04d}.png')
            plt.savefig(outfig, transparent=transparent)


if makefig and savefig and not makegif:
    inum = 1
    orig = savefile
    while os.path.exists(savefile):
        savefile = os.path.splitext(orig)[0] + f'{inum}' + os.path.splitext(orig)[1]
        inum += 1
    plt.savefig(savefile, transparent=transparent)


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


"""
corner = 4

import scipy.ndimage
from astropy.modeling import models, fitting

for corner in np.arange(4)+1:
    fdata = data * 1
    fdata[fdata > 1000] = np.median(fdata)
    fdata[fdata < 50] = 50
    fdata = scipy.ndimage.uniform_filter(fdata, size=201)
    
    p_init = models.Polynomial2D(degree=2)
    fit_p = fitting.LevMarLSQFitter()
    
    xl, yl = fdata.shape
    
    if corner==1:
        xs = lat[:xl//4,:yl//4]
        ys = lon[:xl//4,:yl//4]
        dat = fdata[:xl//4,:yl//4]
    elif corner==2:
        xs = lat[-xl//4:,:yl//4]
        ys = lon[-xl//4:,:yl//4]
        dat = fdata[-xl//4:,:yl//4]
    elif corner==3:
        xs = lat[:xl//4,-yl//4:]
        ys = lon[:xl//4,-yl//4:]
        dat = fdata[:xl//4,-yl//4:]
    elif corner==4:
        xs = lat[-xl//4:,-yl//4:]
        ys = lon[-xl//4:,-yl//4:]
        dat = fdata[-xl//4:,-yl//4:]
    else:
        raise Exception('bad corner')
    
    poly = fit_p(p_init, xs, ys, dat)
    
    
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(lat[1:,1:], lon[1:,1:], fdata, cmap='gray',
                       linewidth=0, antialiased=False, alpha=0.2)
    
    mod = poly(xs, ys)
    
    ax.plot_surface(xs, ys, mod, cmap='viridis',
                       linewidth=0, antialiased=False, alpha=0.2)
    
    dat -= mod
    dat += mod.min()
    
    ax.plot_surface(lat[1:,1:], lon[1:,1:], fdata, cmap='gray',
                       linewidth=0, antialiased=False, alpha=0.2)
    
    ax.text(lat[0,0],lon[0,0],50,'C1',color='red')
    ax.text(lat[-1,0],lon[-1,0],50,'C2',color='red')
    ax.text(lat[0,-1],lon[0,-1],50,'C3',color='red')
    ax.text(lat[-1,-1],lon[-1,-1],50,'C4',color='red')
    plt.title('C{0}'.format(corner))

"""

"""
files = glob(os.path.join(datadir, '*-1-3*fits')) + glob(os.path.join(datadir, '*-1-4*fits'))

latedge, lonedge = [], []
for ii, ifile in enumerate(files):
    with fits.open(ifile) as ff:
        print(ii+1, len(files))
        wcs = WCS(ff[1].header)
        
        data = ff[1].data * 1

        xinds = np.arange(-0.5, data.shape[0]-0.4)
        yinds = np.arange(-0.5, data.shape[1]-0.4)
        mesh = np.meshgrid(xinds, yinds, indexing='ij')
        
        lon, lat = wcs.all_pix2world(mesh[1].flatten(), mesh[0].flatten(), 0)
        lon = lon.reshape(mesh[0].shape)
        lat = lat.reshape(mesh[1].shape)
        lon -= 180.
        
        latedge += [lat[0,0], lat[0, -1], lat[-1, 0], lat[-1, -1]]
        lonedge += [lon[0,0], lon[0, -1], lon[-1, 0], lon[-1, -1]]
"""