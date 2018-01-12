import numpy as np
from numpy.testing import assert_almost_equal
import cubedsphere as cs


def test_csgrid_gmap():

    # read reference data
    ds_FV3 = cs.open_FV3data("./example_notebooks/sample_data/FV3diag_C48/",
                             'atmos_daily')

    # calculate grid online
    grid = cs.csgrid_GMAO(48)

    # fix a single polar point before testing. should not have real effect
    index = (5, 24, 24)  # index of the pole point
    assert grid['lat_b'][index] == -90.0
    assert grid['lon_b'][index] == 35.0
    assert ds_FV3['lon_b'][index] == 350.0
    grid['lon_b'][index] = 350.0  # change to FV3 value

    for varname in ['lon', 'lon_b', 'lat', 'lat_b']:
        assert_almost_equal(grid[varname], ds_FV3[varname], decimal=4)
