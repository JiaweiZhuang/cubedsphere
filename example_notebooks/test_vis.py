import xarray as xr
import cubedsphere
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from matplotlib.colors import LogNorm

plt.ion()

ORTHO = ccrs.Orthographic(-95., 40)


ds = xr.open_dataset(
    "../cubed_sphere/TEST7.geosgcm_prog.20000415_0000z.nc4"
)
data = ds['O3'].isel(time=0).squeeze()
data = data.isel(lev=9)
norm = LogNorm(vmin=3e-8, vmax=1e-7)
# Re-order to (lon, lat, face)
data = data.transpose()

cs = cubedsphere.CSGrid(12)
fig = plt.figure()
ax = fig.add_subplot(111, projection=ORTHO)
coll, ax = cubedsphere.plotting.pcolormesh(
    data, cs, ax=ax,
    # vmin=3e-8, vmax=1e-7
    norm=norm
)
ax.coastlines()

plt.show()
