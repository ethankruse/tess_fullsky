def clean_corner(indata, cornerdir, adjustments, cleanplot=False, sec=None,
                 cam=None, ccd=None, create=False):
    """
    Remove the glow from a corner of the CCD using an empirical model for
    that sector. Also used to help create that empirical model.

    Parameters
    ----------
    indata : np.ndarray
        The array to have its corner glow removed.
    cornerdir : path
        Directory where the empirical corner models are stored.
    adjustments : (np.ndarray, np.ndarray, np.ndarray, np.ndarray)
        Tuple containing the (sector, camera, CCD, manual adjustment level)
        for all the manual adjustments needed for every CCD in every sector to
        make for the cleanest output image.
    cleanplot : bool
        Whether to show a diagnostic plot with the original array, the model,
        and the cleaned array.
    sec : int
        What sector the data came from
    cam : int
        What camera the data came from
    ccd : int
        What CCD the data came from
    create : bool
        Are we creating the empirical model for this sector

    Returns
    -------
    np.ndarray
        If create is False, returns the cleaned version of the input CCD with
        the corner glow removed
    or
    (np.ndarray, np.ndarray, np.ndarray)
        If create is True, returns x coordinates, y coordinates, and smoothed
        corner of this CCD.
    """
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.ndimage

    # make a smoother copy of the input data
    fdata = indata * 1
    # remove spiky stars and bad data
    fdata[fdata > 1000] = np.median(fdata)
    fdata[fdata < 50] = 50
    # smooth things out
    fdata = scipy.ndimage.uniform_filter(fdata, size=201)

    # set up the coordinates
    xinds = np.arange(0.5, indata.shape[0] - 0.4)
    yinds = np.arange(0.5, indata.shape[1] - 0.4)
    lon, lat = np.meshgrid(xinds, yinds, indexing='ij')

    # which corner has the glow
    if ccd == 2 or ccd == 4:
        corner = 3
    elif ccd == 1 or ccd == 3:
        corner = 1
    else:
        raise Exception('Bad ccd in clean_corner')

    xl, yl = fdata.shape
    # pick the right corner and get views of the data.
    # cover all cases even though only corners 1 and 3 are real problems with
    # actual TESS data.
    if corner == 1:
        xs = lat[:xl//4, :yl//4]
        ys = lon[:xl//4, :yl//4]
        dat = fdata[:xl//4, :yl//4]
        rdat = indata[:xl//4, :yl//4]
    elif corner == 2:
        xs = lat[-xl//4:, :yl//4]
        ys = lon[-xl//4:, :yl//4]
        dat = fdata[-xl//4:, :yl//4]
        rdat = indata[-xl//4:, :yl//4]
    elif corner == 3:
        xs = lat[:xl//4, -yl//4:]
        ys = lon[:xl//4, -yl//4:]
        dat = fdata[:xl//4, -yl//4:]
        rdat = indata[:xl//4, -yl//4:]
    elif corner == 4:
        xs = lat[-xl//4:, -yl//4:]
        ys = lon[-xl//4:, -yl//4:]
        dat = fdata[-xl//4:, -yl//4:]
        rdat = indata[-xl//4:, -yl//4:]
    else:
        raise Exception(f'Bad corner value {corner}. Must be 1-4.')

    # if we're creating the model, return the coordinates and smoothed version
    # of this corner to process along with all the other CCDs in this sector.
    if create:
        return xs, ys, dat

    # if this corner already should have a model, load it up
    ctxt = os.path.join(cornerdir, f'sector{sec:02d}.corner.txt')
    fix = np.loadtxt(ctxt)

    if corner == 3:
        fix = fix[:, ::-1]
    elif corner == 1:
        pass
    else:
        raise Exception('Bad corner')

    # unpack these again
    bsec, bcam, bccd, badj = adjustments

    # find our manual adjustment to this corner
    srch = np.where((sec == bsec) & (cam == bcam) & (ccd == bccd))[0]
    if len(srch) > 1:
        raise Exception(f'Multiple adjustments for {sec}, {cam}, {ccd}')

    if len(srch) == 0:
        # adj = 0
        raise Exception(f'No adjustment found for {sec}, {cam}, {ccd}')
    else:
        adj = badj[srch[0]]
    # make the adjustment
    fix *= adj

    # diagnostic plots to make sure it's working if desired
    if cleanplot:
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        # plot the smoothed data
        ax.plot_surface(lat, lon, fdata, cmap='gray',
                        linewidth=0, antialiased=False, alpha=0.3)
        # fix the corner and replot the smoothed but fixed version
        dat -= fix
        ax.plot_surface(lat, lon, fdata, cmap='viridis',
                        linewidth=0, antialiased=False, alpha=0.6)

        # label the corners and the plot
        ax.text(lat[0, 0], lon[0, 0], 50, 'C1', color='red')
        ax.text(lat[-1, 0], lon[-1, 0], 50, 'C2', color='red')
        ax.text(lat[0, -1], lon[0, -1], 50, 'C3', color='red')
        ax.text(lat[-1, -1], lon[-1, -1], 50, 'C4', color='red')
        plt.title(f'Sec {sec}, Cam {cam}, CCD {ccd}, Corner {corner}')

    # remove the trend from the actual (unsmoothed) data
    rdat -= fix

    return indata
