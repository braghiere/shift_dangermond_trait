# Import the NetCDF library for accessing and manipulating NetCDF libraries.
import netCDF4 as nc

# Open the Par8 NetCDF file downloaded previously.
dataset = nc.Dataset('output_clipped.nc')

# Print the high-level metadata.
print(dataset)

# Longitude and latitude.
lons = dataset['longitude'][:]
lats = dataset['latitude'][:]

# Temperature at 4pm (zero-based).
refl = dataset['reflectance'][0,40,:,:]
print(f'Variable shape: {refl.shape}')

import numpy as np
import math

# Scan the extracted data for the minimum and maximum.
min_value = math.floor(np.amin(refl))
max_value = math.ceil(np.amax(refl))
print(f'min value: {min_value}')
print(f'max_value: {max_value}')

import regrid.builder as RegridBuilder

# Enable debug logging.
import logging
logging.basicConfig(level=logging.DEBUG)

regridder = RegridBuilder.build_from_input_grid(lats, lons, resolution=0.01)
regridded_refl = regridder.regrid(refl)

# Disable debug logging.
logging.basicConfig(level=logging.WARNING)

print(f'regridded_refl.shape: {np.shape(regridded_refl)}')
