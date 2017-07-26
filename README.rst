CubedSphere: Cubed-Sphere data processing with xarray 
======================================

**This package is still under initial development, not ready for public use.**

**CubedSphere** is a xarray_-based package for processing the 
`Cubed-Sphere <http://acmg.seas.harvard.edu/geos/cubed_sphere.html>`_ data from  
`GFDL-FV3 <https://www.gfdl.noaa.gov/fv3/>`_ and  
`models using FV3 <https://www.gfdl.noaa.gov/fv3/fv3-applications/>`_ such as  
`GEOS-Chem HP <http://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-Chem_HP>`_,  
`GEOS5 <https://gmao.gsfc.nasa.gov/GEOS/>`_,  
and `CESM <http://www.cesm.ucar.edu>`_.

Installation
------------

Install dependencies via `conda <https://www.continuum.io/downloads>`_::

    $ conda install xarray dask 

Install this package::

    $ git clone http://github.com/jiaweizhuang/cubedsphere
    $ cd cubedsphere
    $ pip install .

For developers, install by symlink to automatically track source file changes::

    $ pip install -e .

License
---------------
This code is freely available under the MIT license.
See the accompanying LICENSE file.

.. _xarray: http://xarray.pydata.org
