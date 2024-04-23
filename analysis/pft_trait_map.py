import xarray as xr
import numpy as np

# Load chl and PFT datasets
chl_ds = xr.open_dataset('chl_aviris_dangermond_time_00_reg.nc')
lai_ds = xr.open_dataset('lai_aviris_dangermond_time_00_reg.nc')
lma_ds = xr.open_dataset('lma_aviris_dangermond_time_00_reg.nc')
pft_ds = xr.open_dataset('../California_Vegetation_WHRTYPE_Dangermond/output_latlon.nc')

# Extract chl values and PFT data
chl_values = chl_ds['chl'].values
lai_values = lai_ds['lai'].values
lma_values = lma_ds['lma'].values
pft_values = pft_ds['Band1'].values  # Assuming PFT values are stored in Band1 variable

# Define PFT categories (2, 3, 4)
pft_categories = [2, 3, 4]

# Create a mask for PFT categories
pft_mask = np.isin(pft_values, pft_categories)

# Ensure chl_values and pft_mask have the same shape
masked_chl_values = np.ma.masked_where(~pft_mask[:, :], chl_values[0, :, :])
masked_lai_values = np.ma.masked_where(~pft_mask[:, :], lai_values[0, :, :])
masked_lma_values = np.ma.masked_where(~pft_mask[:, :], lma_values[0, :, :])

# Calculate average chl value per PFT
average_chl_per_pft = []
std_chl_per_pft = []
average_lai_per_pft = []
std_lai_per_pft = []
average_lma_per_pft = []
std_lma_per_pft = []
for pft_category in pft_categories:
    #chl
    pft_chl_values = masked_chl_values[pft_values[:, :] == pft_category]
    average_chl = np.nanmean(pft_chl_values)
    std_chl = np.nanstd(pft_chl_values)
    average_chl_per_pft.append(average_chl)
    std_chl_per_pft.append(std_chl)
    #lai
    pft_lai_values = masked_lai_values[pft_values[:, :] == pft_category]
    average_lai = np.nanmean(pft_lai_values)
    std_lai = np.nanstd(pft_lai_values)
    average_lai_per_pft.append(average_lai)
    std_lai_per_pft.append(std_lai)
    #lma
    pft_lma_values = masked_lma_values[pft_values[:, :] == pft_category]
    pft_lma_values = pft_lma_values[(~np.isnan(pft_lma_values)) & (pft_lma_values != 0)]
    average_lma = np.nanmean(pft_lma_values)
    std_lma = np.nanstd(pft_lma_values)
    average_lma_per_pft.append(average_lma)
    std_lma_per_pft.append(std_lma)

# Print average and standard deviation chl values per PFT
for idx, (avg_chl, std_chl) in enumerate(zip(average_chl_per_pft, std_chl_per_pft)):
    print(f'PFT {pft_categories[idx]}: Average chl value: {avg_chl}, Standard Deviation: {std_chl}')

# Print average and standard deviation chl values per PFT
for idx, (avg_lai, std_lai) in enumerate(zip(average_lai_per_pft, std_lai_per_pft)):
    print(f'PFT {pft_categories[idx]}: Average lai value: {avg_lai}, Standard Deviation: {std_lai}')

# Print average and standard deviation chl values per PFT
for idx, (avg_lma, std_lma) in enumerate(zip(average_lma_per_pft, std_lma_per_pft)):
    print(f'PFT {pft_categories[idx]}: Average lma value: {avg_lma}, Standard Deviation: {std_lma}')

# Optionally, save the new chl map with masked values
#chl_ds['chl'][:] = masked_chl_values
#chl_ds.to_netcdf('masked_chl_aviris_dangermond.nc')

# Create a new chl map with average values per PFT
new_chl_values = np.zeros_like(chl_values[0,:,:])
new_std_chl_values = np.zeros_like(chl_values[0,:,:])
new_lai_values = np.zeros_like(lai_values[0,:,:])
new_std_lai_values = np.zeros_like(lai_values[0,:,:])
new_lma_values = np.zeros_like(lma_values[0,:,:])
new_std_lma_values = np.zeros_like(lma_values[0,:,:])

for idx, pft_category in enumerate(pft_categories):
    mask = (pft_values == pft_category)[:, :]
    #print(np.shape(chl_values))
    new_chl_values[mask] = average_chl_per_pft[idx]
    new_std_chl_values[mask] = std_chl_per_pft[idx]
    new_lai_values[mask] = average_lai_per_pft[idx]
    new_std_lai_values[mask] = std_lai_per_pft[idx]
    new_lma_values[mask] = average_lma_per_pft[idx]
    new_std_lma_values[mask] = std_lma_per_pft[idx]

# Create new xarray DataArrays for chl and std
new_chl_da = xr.DataArray(new_chl_values[np.newaxis,:,:], dims=('time', 'lat', 'lon'), coords={'time': chl_ds['time'], 'lat': chl_ds['lat'], 'lon': chl_ds['lon']})
new_std_chl_da = xr.DataArray(new_std_chl_values[np.newaxis,:,:], dims=('time', 'lat', 'lon'), coords={'time': chl_ds['time'], 'lat': chl_ds['lat'], 'lon': chl_ds['lon']})

# Create a new xarray Dataset with chl and std DataArrays
new_chl_data = {'chl': new_chl_da, 'std': new_std_chl_da}
new_chl_ds = xr.Dataset(new_chl_data)

# Save the new chl map with masked values and standard deviation to a NetCDF file
new_chl_ds.to_netcdf('masked_chl_aviris_dangermond.nc')

# Create new xarray DataArrays for chl and std
new_lai_da = xr.DataArray(new_lai_values[np.newaxis,:,:], dims=('time', 'lat', 'lon'), coords={'time': lai_ds['time'], 'lat': lai_ds['lat'], 'lon': lai_ds['lon']})
new_std_lai_da = xr.DataArray(new_std_lai_values[np.newaxis,:,:], dims=('time', 'lat', 'lon'), coords={'time': lai_ds['time'], 'lat': lai_ds['lat'], 'lon': lai_ds['lon']})

# Create a new xarray Dataset with chl and std DataArrays
new_lai_data = {'lai': new_lai_da, 'std': new_std_lai_da}
new_lai_ds = xr.Dataset(new_lai_data)

# Save the new chl map with masked values and standard deviation to a NetCDF file
new_lai_ds.to_netcdf('masked_lai_aviris_dangermond.nc')

# Create new xarray DataArrays for chl and std
new_lma_da = xr.DataArray(new_lma_values[np.newaxis,:,:], dims=('time', 'lat', 'lon'), coords={'time': lma_ds['time'], 'lat': lma_ds['lat'], 'lon': lma_ds['lon']})
new_std_lma_da = xr.DataArray(new_std_lma_values[np.newaxis,:,:], dims=('time', 'lat', 'lon'), coords={'time': lma_ds['time'], 'lat': lma_ds['lat'], 'lon': lma_ds['lon']})

# Create a new xarray Dataset with chl and std DataArrays
new_lma_data = {'lma': new_lma_da, 'std': new_std_lma_da}
new_lma_ds = xr.Dataset(new_lma_data)

# Save the new chl map with masked values and standard deviation to a NetCDF file
new_lma_ds.to_netcdf('masked_lma_aviris_dangermond.nc')


