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
        

#From Cheng et al 2014
dsub1368 = dsub.sel(wavelength=slice(1365, 1371)).mean(dim="wavelength")
b1 = dsub1368.reflectance
dsub1722 = dsub.sel(wavelength=slice(1719, 1725)).mean(dim="wavelength")
b2 = dsub1722.reflectance
ndlma = (b1 - b2)/(b1 + b2) 

print(b1,b2)

lma = -46.03 + 1194.54*(ndlma)

#lma = lma.where(lma > -1000.0, np.nan)

lma = lma.where(lma >= 0.0, 0.0)
lma = lma.where(lma < 1000.0, 1000.)

#from g.m-2 to g.cm-2
lma = lma/10000

#ds = dsub.reflectance
ds = lma
        
 
# Create a new DataArray without 'x' and 'y' dimensions
ds = xr.DataArray(ds.values, coords={'time': ds['time'], 'latitude': ds['latitude'], 'longitude': ds['longitude']},
                  dims=['time', 'latitude', 'longitude'])
        
# Reorder dimensions to 'time', 'wavelength', 'latitude', 'longitude'
ds = ds.transpose('time', 'latitude', 'longitude')
        
#print(ds)
        
ds.name = 'lma'
ds.attrs['units'] = 'g.cm-2'

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
ds.to_netcdf(f'lma_aviris_dangermond_time_00.nc')
                
#print(ds)

                
#print(ds)

