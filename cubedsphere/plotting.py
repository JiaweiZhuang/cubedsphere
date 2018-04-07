import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from matplotlib.path import Path

import cartopy.crs as ccrs

from . util import shift_lons, calc_ϕ_on_great_circle

FACE_ROTS = [None, 0., 0., 90, 180, 180, 90]
FACE_ORDER = [1, 2, 6, 4, 5, 3]

def pcolormesh(data, cs, face_wrap=False, ax=None, **kwargs):
    """ Emulate pcolormesh for a cubed-sphere dataset.

    This method helps to construct the polygons/patches necessary
    for plotting on a cubed-sphere mesh over a geographic transformation
    or coordinate reference system. We have to make a few hard-coded
    assumptions about the orientation and layout of the data we're
    passed, but this still lets us construct the precise mesh
    on which we will plot our dataset.

    Parameters
    ----------
    data : xarray.DataArray
        A DataArray with dimensions ordered (lon, lat, 6). The 3rd dimension
        should correspond to faces on the cubed-sphere mesh, in an
        order consistent with the hard-coded metadata parameters
        `FACE_ROTS` (how to rotate each face to orient it correctly) and
        `FACE_ORDER` (the tiling of the faces)
    cs : CSGrid
        A CSGrid object with pre-computed mesh geometry
    ax : Axis
        An Axis to plot off; if not provided, one will be generated but
        without any geographic information
    **kwargs : dict
        Additional arguments for stlying

    """

    # Read the faces in the correct order and do a 2D rotation to
    # orient correctly. This should be inferrable from the metadata, but for
    # now we just hard-code in the logic from the reorderCS.m routine. Note
    # that we flip the placement of the top/bottom patches (they're given out
    # of order??) and we explicitly map rotations to this face ordering.
    new_data = []
    for f in FACE_ORDER:
        _d = data.sel(nf=f).values
        nrot = FACE_ROTS[int(f)] // 90
        _d = np.rot90(_d, nrot)

        new_data.append(_d)
    data = np.stack(new_data, axis=-1)

    FF = None

    # Construct Path objects corresponding to each grid cell in the
    # cubed-sphere mesh.
    all_paths = []
    all_corners = []
    all_data = []
    for face in range(6):  # Loop over faces
        if (FF is not None) and (face != FF):
            continue

        lons, lats = cs.lon_edge[..., face], cs.lat_edge[..., face]
        lons = shift_lons(lons)
        face_data = data[..., face].ravel()

        c = np.stack((lons, lats), axis=-1)
        face_corners = np.concatenate([
            c[0:-1, 0:-1],
            c[1:, 0:-1],
            c[1:, 1:],
            c[0:-1:, 1:],
            c[0:-1, 0:-1]
        ], axis=2)

        face_corners = face_corners.reshape((cs.c * cs.c, 5, 2))
        face_corners = face_corners.tolist()

        new_face_corners = []
        new_face_data = []
        # COUNT = 0

        for corners, datum in zip(face_corners, face_data):
            ll, ul, ur, lr, _ = corners
            if face_wrap and (np.abs(ul[0]) > 160) and (ul[0] * ur[0]) < 0:
                # If we are dealing with cells which cross the date
                # line, then we will artifically split them and
                # construct cells explicitly on either side of it
                λn = -179.9
                λp = 179.9
                ϕu = calc_ϕ_on_great_circle(180., ul, ur)
                ϕl = calc_ϕ_on_great_circle(180., ll, lr)

                new_face_corners.append([ll, ul, [λn, ϕu], [λn, ϕl], ll])
                new_face_corners.append([[λp, ϕu], [λp, ϕl], lr, ur, [λp, ϕu]])
                new_face_data.append(datum)
                new_face_data.append(datum)

                # COUNT += 1
            else:
                new_face_corners.append(corners)
                new_face_data.append(datum)
                continue

        # print(face, COUNT)

        paths = [Path(x) for x in new_face_corners]

        if new_face_corners:
            all_paths.extend(paths)
            all_corners.append(new_face_corners)
            all_data.append(new_face_data)


    all_corners = np.concatenate(all_corners)
    all_data = np.concatenate(all_data)

    # Create PolyCollection and add to
    transform = ccrs.Geodetic()
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        transform = None

    coll = PolyCollection(
        all_corners, lw=0.5, edgecolors='k',
        transform=transform
    )

    coll.set_array(all_data)
    coll.set_alpha(1.0)
    if 'cmap' in kwargs:
        coll.set_cmap(kwargs['cmap'])
    if 'norm' in kwargs:
        norm = kwargs['norm']
        coll.set_norm(norm)
        coll.set_clim(norm.vmin, norm.vmax)
    elif ('vmin' in kwargs) and ('vmax' in kwargs):
        vmin, vmax = kwargs['vmin'], kwargs['vmax']
        coll.set_norm(Normalize(vmin, vmax))
        coll.set_clim(vmin, vmax)

    ax.add_collection(coll, autolim=False)
    if transform is not None:
        ax.autoscale_view()

    return coll, ax


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
        im = ax.pcolormesh(lon_b[i, :, :], lat_b[i, :, :], data[i, :, :],
                           vmin=vmin, vmax=vmax, **kwargs)

    return im


def plotCS_quick(dr, ds, ax, **kwargs):
    """A quick way to plot cubed-sphere data of resolution 6*N*N.

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
    """A wrapper of plotting methods.

    Just very preliminary. See plotCS_quick for now.

    """

    if method == 'quick':
        plotfunc = plotCS_quick
    else:
        raise ValueError("Unrecognized plotting methods!")

    im = plotfunc(*args, **kwargs)

    return im
