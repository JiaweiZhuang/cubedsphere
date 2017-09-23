CubedSphere: Cubed-Sphere data processing with xarray 
=====================================================

**This package is still under initial development, not ready for public use.**

**CubedSphere** is an xarray_-based package for processing the 
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

Developers can install by symlink to automatically track source file changes::

    $ pip install -e .

Users can also install without cloning the source code::

    $ pip install git+https://github.com/JiaweiZhuang/cubedsphere.git

Quick Start
-----------

See the `example notebook <https://github.com/JiaweiZhuang/cubedsphere/blob/master/example_notebooks/basic_design.ipynb>`_.

License
-------
This code is freely available under the MIT license.
See the accompanying LICENSE file.

.. _xarray: http://xarray.pydata.org
