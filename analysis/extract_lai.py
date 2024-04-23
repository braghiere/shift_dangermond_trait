import xarray as xr
from matplotlib import pyplot as plt
import rioxarray as rxr
import os

#dat = xr.open_dataset('/net/fluo/data1/students/renato/aviris_dangermond/aviris_dangermond_time_00/output_clipped.nc')
dat = xr.open_dataset('/net/fluo/data3/data/FluoData1/students/renato/aviris_dangermond/aviris_dangermond_time_01/output_clipped.nc')

import xarray as xr
import numpy as np
import pyproj
import sys

# Define the UTM and WGS84 projections
utm_proj = pyproj.CRS("EPSG:32610")
wgs84_proj = pyproj.CRS("EPSG:4326")

dsub=dat
        
#From landsat 8
dsub_red = dsub.sel(wavelength=slice(636, 673)).mean(dim="wavelength")
red = dsub_red.reflectance
#dsub_green = dsub.sel(wavelength=slice(530, 590)).mean(dim="wavelength")
#green = dsub_green.reflectance
dsub_nir = dsub.sel(wavelength=slice(851, 879)).mean(dim="wavelength")
nir = dsub_nir.reflectance

#SAVI = 1.1*(red-green)/(0.1 + (red-green))
#LAI = -np.log((0.69 - SAVI)/0.59)/0.91
        
SR = nir/red
#From Blinn et al. (2019)
LAI = 0.332915*SR - 0.00212
                

#ds = dsub.reflectance
ds = LAI
        
 
# Create a new DataArray without 'x' and 'y' dimensions
ds = xr.DataArray(ds.values, coords={'time': ds['time'], 'latitude': ds['latitude'], 'longitude': ds['longitude']},
                  dims=['time', 'latitude', 'longitude'])
        
# Reorder dimensions to 'time', 'wavelength', 'latitude', 'longitude'
ds = ds.transpose('time', 'latitude', 'longitude')
        
#print(ds)
        
ds.name = 'lai'

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
#ds.to_netcdf(f'lai_aviris_dangermond_time_00.nc')
ds.to_netcdf(f'lai_aviris_dangermond_time_01.nc')
#print(ds)

                
#print(ds)

