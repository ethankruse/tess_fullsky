def truncate_colormap(origmap, minval=0.0, maxval=1.0, n=-1):
    """
    Create a custom colormap out of a subset of a default matplotlib map.

    Taken from https://stackoverflow.com/a/18926541

    Parameters
    ----------
    origmap : matplotlib.colors.Colormap
        Original colormap to truncate
    minval: float
        Fraction of the way through the map to begin our subset map (0-1).
    maxval: float
        Fraction of the way through the map to end our subset map (0-1).
    n : int
        Number of interpolations in the output map. If -1, use the same number
        as the input map.

    Returns
    -------
    matplotlib.colors.LinearSegmentedColormap
    """
    import matplotlib.colors as mcolors
    import numpy as np
    if n == -1:
        n = origmap.N
    newname = f'trunc({origmap.name},{minval:.2f},{maxval:.2f})'
    values = origmap(np.linspace(minval, maxval, n))
    new_cmap = mcolors.LinearSegmentedColormap.from_list(newname, values)
    return new_cmap
