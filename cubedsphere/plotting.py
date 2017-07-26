import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

def plotCS_quick_raw(lon_b, lat_b, data, lon, ax,
                     masksize=2, vmin=None, vmax=None, **kwargs):
    """A quick way to plot cubed-sphere data of resolution 6*N*N.
    
    This function assumes raw numpy array. 
    For real work it is recommended to use wrapper functions.
    
    Cells near the map boundary are masked to avoid spanning over the entire map.
    Such simple masking only works with the the PlateCarree projection in cartopy.
    
    (TBD) To correctly plot boundary cells or use other projections, use ...

    Parameters
    ----------
    lon_b : numpy array of shape (6,N+1,N+1)
        Longitute of cell boundaries
    
    lan_b : numpy array of shape (6,N+1,N+1)
        Latitude of cell boundaries
        
    data : numpy array of shape (6,N,N)
        The data to be plotted
        
    lon : numpy array of shape (6,N,N)
        Longitute of cell centers, only for masking    
        
    ax : Cartopy GeoAxes
        The axis to be plotted on. 
        This quick plot method only works with the PlateCarree projection.
        ax = plt.axes(projection = ccrs.PlateCarree())
        
    masksize : scalar (in longitude), optional
        The size of the map boundary mask. 
        You can use a smaller mask at higher resolution.
    
    Other keyword arguments such as cmap will be passed to plt.pcolormesh()
    
    Returns
    -------
    im : matplotlib QuadMesh object
        It contains the color information for making colorbar  
    """
    
    # Have to use the same range for all tiles!
    # cannot just pass None to plt.pcolormesh()
    if vmin is None: 
        vmin = data.min()
    if vmax is None: 
        vmax = data.max()
        
    # mask cells near boundaries, otherwise they will span the entire map 
    mask = np.abs(lon - 180) < masksize    
    data = np.ma.masked_where(mask, data)
    
    for i in range(6):
        # 6 tiles have the same color configuration so we only return one QuadMesh object
        im = ax.pcolormesh(lon_b[i,:,:], lat_b[i,:,:], data[i,:,:],
                           vmin=vmin, vmax=vmax, **kwargs)
        
    return im


def plotCS_quick(dr, ds, ax, **kwargs):
    """ A quick way to plot cubed-sphere data of resolution 6*N*N.
    
    Wrapping plotCS_quick_raw to work with xarray objects
    
    Parameters
    ----------
    dr : xarray DataArray
        The dimensions must be (tile,y,x)
    
    ds : xarray DataSet (the parent DataSet of dr)
        Must contain lon_b, lat_b, lon as coordinate variables.
                
        TBD: provide a cleaner API so we can pass a single object.
        
    ax : Cartopy GeoAxes
        The axis to be plotted on. 
        This quick plot method only works with the PlateCarree projection.
        ax = plt.axes(projection = ccrs.PlateCarree())
        
    Other keyword arguments such as cmap will be passed to plotCS_quick_raw(),
    then to plt.pcolormesh()
    
    Returns
    -------
    im : matplotlib QuadMesh object
        It contains the color information for making colorbar  
    """
    
    # must convert xarray objects to raw numpy arrays
    # otherwise numpy masking functions won't work
    
    im = plotCS_quick_raw(ds['lon_b'].values, ds['lat_b'].values,
                          dr.values, ds['lon'].values,
                          ax, **kwargs)
    return im

def plotCS(*args, method='quick', **kwargs):
    """ A wrapper of plotting methods. Just very preliminary.
        See plotCS_quick for now.
    """
    
    if method == 'quick':
        plotfunc = plotCS_quick
    else:
        raise ValueError("Unrecognized plotting methods!")
    
    im = plotfunc(*args, **kwargs)
    
    return im

