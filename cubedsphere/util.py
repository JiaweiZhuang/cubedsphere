import numpy as np
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
