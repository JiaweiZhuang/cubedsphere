from setuptools import setup

setup(name='cubedsphere',
      version='0.0.1',
      description='Cubed-Sphere data processing with xarray',
      url='http://github.com/jiaweizhuang/cubedsphere',
      author='Jiawei Zhuang',
      author_email='jiaweizhuang@g.harvard.edu',
      license='MIT',
      packages=['cubedsphere'],
      install_requires=['xarray','dask','numpy','matplotlib','cartopy']
      )
