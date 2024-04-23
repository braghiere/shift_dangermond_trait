import xarray as xr
import pyproj
import numpy as np
from datetime import datetime

# Define the UTM and WGS84 projections
utm_proj = pyproj.CRS("EPSG:32610")
wgs84_proj = pyproj.CRS("EPSG:4326")

# Define the correct dates for each time folder
dates = [
    '2022-02-24T00:00:00.000000', '2022-02-28T00:00:00.000000',
    '2022-03-08T00:00:00.000000', '2022-03-16T00:00:00.000000',
    '2022-03-22T00:00:00.000000', '2022-04-05T00:00:00.000000',
    '2022-04-12T00:00:00.000000', '2022-04-20T00:00:00.000000',
    '2022-04-29T00:00:00.000000', '2022-05-03T00:00:00.000000',
    '2022-05-11T00:00:00.000000', '2022-05-17T00:00:00.000000',
    '2022-05-29T00:00:00.000000'
]

# Loop over folders from aviris_dangermond_time_00 to aviris_dangermond_time_12
for folder_num in range(13):
    folder = f'/net/fluo/data3/data/FluoData1/students/renato/aviris_dangermond/aviris_dangermond_time_{folder_num:02d}'
    print(f'Processing data in folder: {folder}...')
    
    # Load hyperspectral data from the clipped netCDF file
    dat = xr.open_dataset(f'{folder}/output_clipped.nc')
    
    # Calculate LAI
    b1368 = dat.sel(wavelength=slice(1365, 1371)).mean(dim="wavelength").reflectance
    b1722 = dat.sel(wavelength=slice(1719, 1725)).mean(dim="wavelength").reflectance
    ndlma = (b1368 - b1722) / (b1368 + b1722)
    lma = -46.03 + 1194.54 * ndlma
    #lma = np.clip(lma, 0, 1000)  # Clip values between 0 and 1000
    lma = np.clip(lma, 25, 200)  # Clip values between 25 and 200
    
    # Calculate CHL
    green_band = dat.sel(wavelength=slice(540, 560)).mean(dim="wavelength").reflectance
    nir_band = dat.sel(wavelength=slice(760, 800)).mean(dim="wavelength").reflectance
    chl = (1 / green_band - 1 / nir_band) * nir_band
    chl = chl * 0.089348898  # Convert to ug/cm^2
    # Replace negative values with NaNs
    chl = chl.where(chl > 0)
    chl = np.clip(chl, 1e-3, float(np.max(chl)))  # Clip values between 0 and max_chl
    
    # Calculate LAI from hyperspectral data
    red = dat.sel(wavelength=slice(636, 673)).mean(dim="wavelength").reflectance
    nir = dat.sel(wavelength=slice(851, 879)).mean(dim="wavelength").reflectance
    SR = nir / red
    lai = 0.332915 * SR - 0.00212
    # Replace negative values with NaNs
    lai = lai.where(lai > 0)
    lai = np.clip(lai, 1e-3, 20.)  # Clip values between 0 and max_lai
    
    # Prepare the DataArrays for writing to netCDF
    lma_ds = xr.DataArray(
        lma,
        dims=('time', 'latitude', 'longitude'),
        coords={'time': dat['time'], 'latitude': dat['latitude'], 'longitude': dat['longitude']},
    )
    lma_ds['latitude'].attrs['standard_name'] = 'latitude'
    lma_ds['longitude'].attrs['standard_name'] = 'longitude'
    lma_ds['time'] = np.datetime64(datetime.strptime(dates[folder_num], '%Y-%m-%dT%H:%M:%S.%f'))
    lma_ds.name = 'lma'
    
    chl_ds = xr.DataArray(
        chl,
        dims=('time', 'latitude', 'longitude'),
        coords={'time': dat['time'], 'latitude': dat['latitude'], 'longitude': dat['longitude']},
    )
    chl_ds['latitude'].attrs['standard_name'] = 'latitude'
    chl_ds['longitude'].attrs['standard_name'] = 'longitude'
    chl_ds['time'] = np.datetime64(datetime.strptime(dates[folder_num], '%Y-%m-%dT%H:%M:%S.%f'))
    chl_ds.name = 'chl'
    
    lai_ds = xr.DataArray(
        lai,
        dims=('time', 'latitude', 'longitude'),
        coords={'time': dat['time'], 'latitude': dat['latitude'], 'longitude': dat['longitude']},
    )
    lai_ds['latitude'].attrs['standard_name'] = 'latitude'
    lai_ds['longitude'].attrs['standard_name'] = 'longitude'
    lai_ds['time'] = np.datetime64(datetime.strptime(dates[folder_num], '%Y-%m-%dT%H:%M:%S.%f'))
    lai_ds.name = 'lai'
    
    # Write to netCDF files
    lma_ds.to_netcdf(f'lma_aviris_dangermond_time_{folder_num:02d}_v2.nc')
    chl_ds.to_netcdf(f'chl_aviris_dangermond_time_{folder_num:02d}.nc')
    lai_ds.to_netcdf(f'lai_aviris_dangermond_time_{folder_num:02d}.nc')
    
    print(f'LMA, CHL, and LAI data processed and saved for {dates[folder_num]} in folder: {folder}')

