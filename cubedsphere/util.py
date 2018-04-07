import numpy as np
#from scipy.optimize import minimize_scalar
#from pysal.cg.sphere import geointerpolate, arcdist

import warnings


def _bin_phydim(dr, dim_name, bin_bound, bin_center, stacked_dim_name):
    """Binning over a physical dimension, weighted by area. Internal function
    for taking zonal or meridonal average.

    Parameters
    ----------
    dr : xarray DataArray
        Must contain (tile,y,x) as dimensions.
        Can have arbitrary additional dimensions.

    dim_name : str
        support 'lat' or 'lon'

    bin_bound : 1D array (in degree)
        The boundaries of bins in longitude or latitude

    bin_center : 1D array (in degree)
        The center of bins in longitude or latitude

    stacked_dim_name : str
        Will be 'stacked_tile_y_x' if dimension names are (tile,y,x)

    Returns
    -------
    dr_binned : xarray DataArray

    """

    # use multi-dimensional groupby, with 1D bins
    group = dr.groupby_bins(dim_name, bin_bound, labels=bin_center)

    # select the algorithm depending on the existence of 'area'
    if 'area' in dr:
        # take weighted average

        # a MapReduce algorithm that minimizes the "Reduction" operations
        # that's the fastest way I can find for taking weighted average
        # It is 10x faster than the algorithm on xarray online docs:
        # http://xarray.pydata.org/en/stable/examples/monthly-means.html
        def mean_weighted_by_area(dr):
            weights = dr['area']/dr['area'].sum()
            dr_mean = (dr*weights).sum(dim=stacked_dim_name)
            return dr_mean

        dr_binned = group.apply(mean_weighted_by_area)

    else:
        # simply take unweighted average
        warnings.warn("Use uniform weights because the input DataArray "
                      "does not have area as a coordinate variable. "
                      "The result could be inaccurate because "
                      "Cubed-sphere box size varies."
                      )
        dr_binned = group.mean(dim=stacked_dim_name)

    return dr_binned


def zonal_mean(dr, dlat=4):
    """Take the zonal mean of cubed-sphere data.

    Parameters
    ----------
    dr : xarray DataArray
        Must contain (tile,y,x) as dimensions.
        Can have arbitrary additional dimensions.

    dlat : scalar (in latitude), optional
        The binning size in latitude

    Returns
    -------
    zonal_mean : xarray DataArray

    """

    dim_name = 'lat'
    bin_bound = np.arange(-90, 91, dlat)
    bin_center = np.arange(-89, 90, dlat)
    stacked_dim_name = 'stacked_tile_y_x'

    zonal_mean = _bin_phydim(dr, dim_name, bin_bound, bin_center, stacked_dim_name)

    return zonal_mean


def meridional_mean(dr, dlon=5):
    """Take the meridional mean of cubed-sphere data.

    Parameters
    ----------
    dr : xarray DataArray
        Must contain (tile,y,x) as dimensions.
        Can have arbitrary additional dimensions.

    dlon : scalar (in longitude), optional
        The binning size in latitude

    Returns
    -------
    meri_mean : xarray DataArray

    """

    dim_name = 'lon'
    bin_bound = np.arange(0, 361, dlon)
    bin_center = np.arange(1, 360, dlon)
    stacked_dim_name = 'stacked_tile_y_x'

    meri_mean = _bin_phydim(dr, dim_name, bin_bound, bin_center, stacked_dim_name)

    return meri_mean


def gaussian_bell(xs, ys, xc=0., yc=0., xsigma=1., ysigma=1.):
    """ Compute a 2D Gaussian with asymmetric standard deviations and
    arbitrary center.

    .. math::

        Z = \exp{\left[\frac{(x - x_c)^2}{2\sigma_x} + \frac{(y - y_c)^2}{2\sigma_y}\right]}

    Parameters
    ----------
    {x,y}s : array-like of floats
        x- and y-coordinates where the function should be calculated. Can be
        arbitrary shape as long as they both match.
    {x,y}c : float
        coordinates corresponding to center of bell.
    {x,y}sigma : float
        width/standard deviation (σ) of distribution in each coordinate direction.

    Returns
    -------
    Z evaluated at the given coordinates.

    """
    expon = ((xs - xc)**2)/2./xsigma + ((ys - yc)**2)/2./ysigma
    return np.exp(-expon)

def multi_wave(lons, lats, nx=2, ny=1):
    """ Compute an arbitrary zonally/meridionally varying wave.

    .. math::

        Z = \cos{\frac{n_x \lambda}{T_\lambda}} + 2\frac{\phi - \bar{\phi}}{\mathrm{std}(\phi)}

    Parameters
    ----------
    lons, lats : array-like of floats
        Longitude/latitude coordinate at which to evaluate wave equation
    nx, ny : int
        Wavenumber in zonal and meridional direction

    """
    Tx = 360. / 2. * np.pi
    Ty = 180. / 2. * np.pi
    # return np.sin(nx*lons/Tx + np.cos(lats/Ty)) #+ np.cos(ny*lats/Ty)
    return np.cos(nx*lons/Tx) + 2*(lats - np.mean(lats))/lats.std()


def calc_ϕ_on_great_circle(λp, pt1, pt2, verbose=False):
    """ Calculate latitude corresponding to a given longitude on
    a great circle that passes through two given points.

    Parameters
    ----------
    λp : float
        Longitude to inspect great circle, in degrees
    pt1, pt2 : tuples of (float, float)
        Longitude/latitude of points through which to draw a great
        circle. Coordinates should be passed in degrees
    verbose : logical (default = False)
        Print the result before returning

    Returns
    -------
    latitude corresponding to λp along the desired great circle, in
    degrees

    """
    # distance = arcdist(pt0, pt1, 1.)
    diff = lambda x: np.abs(λp - geointerpolate(pt1, pt2, x)[0])
    result = minimize_scalar(diff, bounds=[0., 1.], method='bounded')
    if verbose:
        print(result)
    λ, ϕ = geointerpolate(pt1, pt2, result.x)
    return ϕ


def shift_lons(lons):
    """ Shift longitudes from (0, 360) to (-180, 180) """
    new_lons = np.empty_like(lons)
    mask = lons > 180
    new_lons[mask] = -(360. - lons[mask])
    new_lons[~mask] = lons[~mask]
    lons = new_lons.copy()
    return lons
