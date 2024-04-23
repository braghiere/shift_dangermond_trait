import xarray as xr
import pyproj

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
    
    # Load LAI data from the clipped netCDF file
    lai_data = xr.open_dataset(f'{folder}/output_clipped.nc')
    
    # Calculate LAI from hyperspectral data
    red = lai_data.sel(wavelength=slice(636, 673)).mean(dim="wavelength").reflectance
    nir = lai_data.sel(wavelength=slice(851, 879)).mean(dim="wavelength").reflectance
    SR = nir / red
    LAI = 0.332915 * SR - 0.00212
    
    # Prepare the DataArray for writing to netCDF
    lai_ds = LAI.transpose('time', 'latitude', 'longitude')
    lai_ds['time'] = lai_ds['time'].astype('int64') / 1e9  # Convert to POSIX timestamp
    lai_ds['time'].attrs['units'] = 'seconds since ' + dates[folder_num]
    lai_ds['time'].attrs['calendar'] = 'gregorian'
    lai_ds.name = 'lai'
    
    # Update coordinate attributes
    lai_ds['latitude'].attrs['standard_name'] = 'latitude'
    lai_ds['longitude'].attrs['standard_name'] = 'longitude'
    
    # Add CRS information
    crs = pyproj.CRS.from_epsg(32610)
    lai_ds.attrs['crs'] = crs
    # Convert crs attribute to string
    lai_ds.attrs['crs'] = str(lai_ds.attrs['crs'])
    
    # Write to netCDF file
    lai_ds.to_netcdf(f'lai_aviris_dangermond_time_{folder_num:02d}.nc')
    
    print(f'LAI data processed and saved for {dates[folder_num]} in folder: {folder}')

