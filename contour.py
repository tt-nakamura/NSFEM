import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.path import Path

N = 100 # number of grid points in x direction

def contour(x, y, z, *a, **k):
    """ plot contour from unstructured data
    x : 1d-array, x coordinates of data points
    y : 1d-array, y coordinates of data points
    z : 1d-array, z coordinates to be plotted
    plotting region is a rectangle given by
      (min(x), max(x), min(y), max(y))
    a : arguments passed to plt.contour
    k : keyword arguments passed to plt.contour
    if k has a key 'confine', then
      k['confine'] : 2d-array, shape(n,2)
        vertices of polygon, outside of which is
        excluded from plotting region
      k['confine'][:,0] = x coordinates of vertices
      k['confine'][:,1] = y corrdinates of vertices
    if k has a key 'exclude', then
      same as 'confine' except that
      inside of polygon is excluded
    """
    x1,x2 = np.min(x), np.max(x)
    y1,y2 = np.min(y), np.max(y)
    X = np.linspace(x1, x2, N)
    Y = np.linspace(y1, y2, N*(y2-y1)/(x2-x1))
    X,Y = np.meshgrid(X,Y)
    XY = np.c_[X.flat, Y.flat]
    Z = griddata(np.c_[x,y], z, XY)
    if 'confine' in k:
        p = Path(k.pop('confine'))
        D = p.contains_points(XY)
        Z[~D] = np.nan
    if 'exclude' in k:
        p = Path(k.pop('exclude'))
        D = p.contains_points(XY)
        Z[D] = np.nan

    Z = Z.reshape(X.shape)
    plt.contour(X, Y, Z, *a, **k)
