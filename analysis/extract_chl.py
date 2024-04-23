import xarray as xr
from matplotlib import pyplot as plt
import rioxarray as rxr
import os

dat = xr.open_dataset('/net/fluo/data1/students/renato/aviris_dangermond/aviris_dangermond_time_00/output_clipped.nc')

import xarray as xr
import numpy as np
import pyproj
import sys

# Define the UTM and WGS84 projections
utm_proj = pyproj.CRS("EPSG:32610")
wgs84_proj = pyproj.CRS("EPSG:4326")

dsub=dat

#From Schneider et al. 
dsub540 = dsub.sel(wavelength=slice(540, 560)).mean(dim="wavelength")
green = dsub540.reflectance
dsub760 = dsub.sel(wavelength=slice(760, 800)).mean(dim="wavelength")
nir = dsub760.reflectance
chl = (1/green - 1/nir) * nir

#From Gitelson et al 2003 (Fig8)
dsub750 = dsub.sel(wavelength=slice(750, 800)).mean(dim="wavelength")
r750 = dsub750.reflectance
dsub695 = dsub.sel(wavelength=slice(695, 740)).mean(dim="wavelength")
r695 = dsub695.reflectance

gitelson2 = (r750/r695) - 1

#In umol.m-2
chl = (800./6.)*gitelson2      
#chl = (1)*gitelson2

#In ug.cm-2
chl = chl*0.089348898

chl = chl.where(chl > 0.0, np.nan)             

#ds = dsub.reflectance
ds = chl
        
 
# Create a new DataArray without 'x' and 'y' dimensions
ds = xr.DataArray(ds.values, coords={'time': ds['time'], 'latitude': ds['latitude'], 'longitude': ds['longitude']},
                  dims=['time', 'latitude', 'longitude'])
        
# Reorder dimensions to 'time', 'wavelength', 'latitude', 'longitude'
ds = ds.transpose('time', 'latitude', 'longitude')
        
#print(ds)
        
ds.name = 'chl'
#ds.attrs['units'] = 'micromol.m-2'
ds.attrs['units'] = 'ug.cm-2'

# Update coordinate attributes
ds['latitude'].attrs['standard_name'] = 'latitude'
ds['longitude'].attrs['standard_name'] = 'longitude'
        

crs = pyproj.CRS.from_epsg(32610)
ds.attrs['crs'] = crs

# Convert crs attribute to string
ds.attrs['crs'] = str(ds.attrs['crs'])

ds['time'] = ds['time'].astype('int64') / 1e9  # Convert to POSIX timestamp
ds['time'].attrs['units'] = 'seconds since 2022-02-24 00:00:00'
ds['time'].attrs['calendar'] = 'gregorian'

# Write to netCDF file
ds.to_netcdf(f'chl_aviris_dangermond_time_00.nc')
                
#print(ds)

                
#print(ds)

