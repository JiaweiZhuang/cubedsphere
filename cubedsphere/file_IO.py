import xarray as xr

def open_FV3data(maindir, prefix):
    """Open output data from the GFDL-FV3 model
    
    All FV3 output data have the name xxx.tile1.nc ~ xxx.tile6.nc, for example
    - atmos_static.tile1.nc ~ atmos_static.tile6.nc
    - atmos_daily.tile1.nc ~ atmos_daily.tile6.nc
    "atmos_static" should always exist and contain the grid information
    
    Parameters
    ----------
    maindir : str
        The directory containing FV3 output files.
    
    prefix : str
        The prefix of FV3 output tile files. For example, "atmos_daily"
    
    Returns
    -------
    ds : xarray DataSet
    """

    # ds will be automatically chunked by tiles. 
    # Not sure if that's the best chunking option.
    ds = xr.open_mfdataset(maindir+prefix+'.tile*.nc',
                           concat_dim='tile', decode_times=False)
    
    grid = xr.open_mfdataset(maindir+'atmos_static.tile*.nc',
                             concat_dim='tile', decode_times=False)
    
    # pass the grid information to the data object, 
    # and use clearer variable names
    ds['lon'] = grid['grid_lont']
    ds['lat'] = grid['grid_latt']
    ds['lon_b'] = grid['grid_lon']
    ds['lat_b'] = grid['grid_lat']
    ds['area'] = grid['area']
    
    # use clearer dimension names
    ds.rename(dict(grid_xt='x', grid_yt='y', 
                   grid_x='x_b', grid_y='y_b'), inplace=True)

    # from data variable to coordinate variable
    ds.set_coords(['lon', 'lat', 'lon_b', 'lat_b','area'], inplace=True)

    return ds

def open_GCHPdata(filename):
    # There's no grid information in GCHP data 
    # so we need to calculate it on-the-fly
    raise NotImplementedError