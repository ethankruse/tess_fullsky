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
cornerdir = os.path.join(os.path.split(__file__)[0], 'corners')

#cS12 2c1 c3c2
#files = glob(os.path.join(datadir, '*s0012-2-1*fits'))
#files = glob(os.path.join(datadir, '*s0001-*fits'))
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

makecorner = False
cornersec = 13

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
transparent = False

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

if makecorner:
    files = glob(os.path.join(datadir, f'*s00{cornersec:02d}-*fits'))
    makefig = False
    savefig = False
    makegif = False
    doclean = True
    test = False

if makegif:
    figdir = os.path.join(figdir, 'gif2')
    prev = glob(os.path.join(figdir, '*png'))
    for iprev in prev:
        os.remove(iprev)
    
if not os.path.exists(figdir):
    os.makedirs(figdir, exist_ok=True)

if test:
    # XXX: for testing
    # files = [files[0]]
    files = files[::3]
    pass


#bsec, bcam, bccd, bcor = np.loadtxt(badcornerfile, unpack=True, ndmin=2, 
#                                    delimiter=',', dtype=int)



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


def clean(data, cleanplot=False, ccd=None, sec=None, cam=None, makecorner=False):
    import scipy.ndimage
    
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
    
    if ccd == 2 or ccd == 4:
        corner = 3
    elif ccd == 1 or ccd == 3:
        corner = 1
    else:
        raise Exception('bad ccd in clean')
        
    # fix the desired corners
    xl, yl = fdata.shape
    
    # pick the right corner and get views of the data
    if corner==1:
        xs = lat[:xl//4,:yl//4]
        ys = lon[:xl//4,:yl//4]
        dat = fdata[:xl//4,:yl//4]
        rdat = data[:xl//4,:yl//4]        
    elif corner==2:
        xs = lat[-xl//4:,:yl//4]
        ys = lon[-xl//4:,:yl//4]
        dat = fdata[-xl//4:,:yl//4]
        rdat = data[-xl//4:,:yl//4]
    elif corner==3:
        xs = lat[:xl//4,-yl//4:]
        ys = lon[:xl//4,-yl//4:]
        dat = fdata[:xl//4,-yl//4:]
        rdat = data[:xl//4,-yl//4:]
    elif corner==4:
        xs = lat[-xl//4:,-yl//4:]
        ys = lon[-xl//4:,-yl//4:]
        dat = fdata[-xl//4:,-yl//4:]
        rdat = data[-xl//4:,-yl//4:]
    else:
        raise Exception('bad corner')
    
    if makecorner:
        return xs, ys, dat
    
    ctxt = os.path.join(cornerdir, f'sector{sec:02d}.corner.txt')
    fix = np.loadtxt(ctxt)
    
    if corner == 3:
        fix = fix[:, ::-1]
    elif corner == 1:
        pass
    else:
        raise Exception('bad corner')
    
    # diagnostic plots to make sure it's working
    if cleanplot:
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot_surface(lat, lon, fdata, cmap='gray',
                           linewidth=0, antialiased=False, alpha=0.3)

        dat -= fix
        
        ax.plot_surface(lat, lon, fdata, cmap='viridis',
                           linewidth=0, antialiased=False, alpha=0.6)
        
        ax.text(lat[0,0],lon[0,0],50,'C1',color='red')
        ax.text(lat[-1,0],lon[-1,0],50,'C2',color='red')
        ax.text(lat[0,-1],lon[0,-1],50,'C3',color='red')
        ax.text(lat[-1,-1],lon[-1,-1],50,'C4',color='red')
        plt.title(f'Sec {sec}, Cam {cam}, CCD {ccd}, Corner {corner}')
    
    # remove the trend from the actual data
    rdat -= fix

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


xxs, yys, dats, ccds = [], [], [], []

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
        # this should always be true
        data = data[:2048, 44:2092]
        lon = lon[:2049, 44:2093]
        lat = lat[:2049, 44:2093]
        
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
            #tocor = np.where((bsec==isec) & (bcam==icam) & (bccd==iccd))[0]
            #print(f'Correcting corners {bcor[tocor]}')
            if makecorner:
                xs, ys, dat = clean(data, cleanplot=cleanplot, makecorner=True,
                                    ccd=iccd, sec=isec, cam=icam)
                xxs.append(xs)
                yys.append(ys)
                dats.append(dat)
                ccds.append(iccd)
            else:
                data = clean(data, cleanplot=cleanplot, ccd=iccd, sec=isec, cam=icam)
            
        if makefig:
            if ii == 0 or test:
                if highres:
                    fig = plt.figure(figsize=(100,100))
                else:
                    fig = plt.figure()
                ax = plt.axes([0.01, 0.01, 0.98, 0.98], projection=tr)
                if highres:
                    ax.outline_patch.set_linewidth(16)
                if transparent:
                    #ax.outline_patch.set_alpha(0)
                    ax.background_patch.set_alpha(0)
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
                outfig = os.path.join(figdir, f'transp_img{(ii+1)//16:04d}.png')
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
        #xs = xxs[ii]
        #ys = yys[ii]
        dat = dats[ii] * 1
        ccd = ccds[ii]
        
        if ccd == 2 or ccd == 4:
            corner = 3
        elif ccd == 1 or ccd == 3:
            corner = 1
        else:
            raise Exception()
            
        if corner == 3:
            dat = dat[:, ::-1]
        
        dx, dy = dat.shape
        
        imin = np.median(dat[3*dx//4:, 3*dy//4])
        ax.plot_surface(xs, ys, dat - imin, cmap='gray',
                        linewidth=0, antialiased=False, alpha=0.2)
        
        stack.append(dat - imin)
        
    avg = np.median(np.dstack(stack), axis=-1)
    ax.plot_surface(xs, ys, avg, cmap='viridis',
                        linewidth=0, antialiased=False, alpha=0.4)
    plt.title(f'Sec {cornersec} Corners')
    
    ctxt = os.path.join(cornerdir, f'sector{cornersec:02d}.corner.txt')
    cfig = os.path.join(cornerdir, f'sector{cornersec:02d}.corner.png')
    np.savetxt(ctxt, avg)
    plt.savefig(cfig)
    