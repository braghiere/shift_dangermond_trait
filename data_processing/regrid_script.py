import xarray as xr
import numpy as np
from scipy.interpolate import griddata

# Function to regrid dataset
def regrid(ds_source, ds_target_grid):
    # Extract source grid
    lon_source, lat_source = np.meshgrid(ds_source['lon'], ds_source['lat'])
    source_grid = np.array([lon_source.flatten(), lat_source.flatten()]).T

    # Extract target grid
    lon_target, lat_target = np.meshgrid(ds_target_grid['lon'], ds_target_grid['lat'])
    target_grid = np.array([lon_target.flatten(), lat_target.flatten()]).T

    # Create a new dataset to hold the regridded data
    ds_regridded = xr.Dataset()

    # Loop over all data variables in the source dataset
    for var in ds_source.data_vars:
        if {'lat', 'lon', 'time'}.issubset(ds_source[var].dims):
            print(f"Regridding {var}")
            # Initialize an empty array to store regridded data
            regridded_data = np.empty((ds_source['time'].size, ds_target_grid['lat'].size, ds_target_grid['lon'].size))

            # Interpolate for each time step
            for t in range(ds_source['time'].size):
                # Interpolate using griddata
                interpolated_data = griddata(source_grid, ds_source[var][t].values.flatten(), target_grid, method='linear')
                # Reshape and assign to the corresponding time step
                regridded_data[t] = interpolated_data.reshape(len(ds_target_grid['lat']), len(ds_target_grid['lon']))

            # Assign to new dataset
            ds_regridded[var] = (('time', 'lat', 'lon'), regridded_data)

    # Assign coordinates to the new dataset
    ds_regridded = ds_regridded.assign_coords({'time': ds_source['time'], 'lat': ds_target_grid['lat'], 'lon': ds_target_grid['lon']})

    return ds_regridded

# Load the source dataset
ds_source = xr.open_dataset('TROPOMI_SIF740nm-v1.001deg_regrid_Dangermond_tll_spatial_ref.nc')

# Load or define the target grid dataset
ds_target_grid = xr.open_dataset('TROPOMI_SIF740nm-v1.001deg_regrid_Dangermond_tll_clipped_458_492.nc')

# Perform regridding
ds_regridded = regrid(ds_source, ds_target_grid)


# Save the regridded dataset to a new NetCDF file
ds_regridded.to_netcdf('TROPOMI_SIF740nm-v1.001deg_regrid_Dangermond_tll_458_492.nc')

print("Regridding complete. Output saved to 'regridded_data.nc'")

