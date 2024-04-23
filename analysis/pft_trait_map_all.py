import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# Define the range of times from 00 to 12
times = ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

# Given dates
dates = ["2022-02-24T00:00:00.000000", "2022-02-28T00:00:00.000000", "2022-03-08T00:00:00.000000",
         "2022-03-16T00:00:00.000000", "2022-03-22T00:00:00.000000", "2022-04-05T00:00:00.000000",
         "2022-04-12T00:00:00.000000", "2022-04-20T00:00:00.000000", "2022-04-29T00:00:00.000000",
         "2022-05-03T00:00:00.000000", "2022-05-11T00:00:00.000000", "2022-05-17T00:00:00.000000",
         "2022-05-29T00:00:00.000000"]

# Extract dates without times
dates_without_times = [date.split('T')[0] for date in dates]

# Loop over times
for time in times:

    # Informative message: Loading datasets for the current time
    print(f'Processing data for time: {time}')
    
    # Load chl and PFT datasets for the current time
    chl_ds = xr.open_dataset(f'chl_aviris_dangermond_time_{time}_reg.nc')
    lai_ds = xr.open_dataset(f'lai_aviris_dangermond_time_{time}_reg.nc')
    lma_ds = xr.open_dataset(f'lma_aviris_dangermond_time_{time}_reg.nc')
    pft_ds = xr.open_dataset(f'../California_Vegetation_WHRTYPE_Dangermond/output_latlon.nc')  # Update the PFT dataset filename
    
    # Informative message: Extracting data from the loaded datasets
    print('Extracting data from the loaded datasets...')
    
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
        # chl
        pft_chl_values = masked_chl_values[pft_values[:, :] == pft_category]
        average_chl = np.nanmean(pft_chl_values)
        std_chl = np.nanstd(pft_chl_values)
        average_chl_per_pft.append(average_chl)
        std_chl_per_pft.append(std_chl)
        
        # lai
        pft_lai_values = masked_lai_values[pft_values[:, :] == pft_category]
        average_lai = np.nanmean(pft_lai_values)
        std_lai = np.nanstd(pft_lai_values)
        average_lai_per_pft.append(average_lai)
        std_lai_per_pft.append(std_lai)
        
        # lma
        pft_lma_values = masked_lma_values[pft_values[:, :] == pft_category]
        pft_lma_values = pft_lma_values[(~np.isnan(pft_lma_values)) & (pft_lma_values != 0)]
        average_lma = np.nanmean(pft_lma_values)
        std_lma = np.nanstd(pft_lma_values)
        average_lma_per_pft.append(average_lma)
        std_lma_per_pft.append(std_lma)
    
    # Print average and standard deviation chl values per PFT
    for idx, (avg_chl, std_chl) in enumerate(zip(average_chl_per_pft, std_chl_per_pft)):
        print(f'PFT {pft_categories[idx]}: Average chl value: {avg_chl}, Standard Deviation: {std_chl}')

    # Print average and standard deviation lai values per PFT
    for idx, (avg_lai, std_lai) in enumerate(zip(average_lai_per_pft, std_lai_per_pft)):
        print(f'PFT {pft_categories[idx]}: Average lai value: {avg_lai}, Standard Deviation: {std_lai}')

    # Print average and standard deviation lma values per PFT
    for idx, (avg_lma, std_lma) in enumerate(zip(average_lma_per_pft, std_lma_per_pft)):
        print(f'PFT {pft_categories[idx]}: Average lma value: {avg_lma}, Standard Deviation: {std_lma}')

    # Optionally, save the new chl map with masked values
    # chl_ds['chl'][:] = masked_chl_values
    # chl_ds.to_netcdf('masked_chl_aviris_dangermond.nc')

    # Create a new chl map with average values per PFT
    new_chl_values = np.zeros_like(chl_values[0,:,:])
    new_std_chl_values = np.zeros_like(chl_values[0,:,:])
    new_lai_values = np.zeros_like(lai_values[0,:,:])
    new_std_lai_values = np.zeros_like(lai_values[0,:,:])
    new_lma_values = np.zeros_like(lma_values[0,:,:])
    new_std_lma_values = np.zeros_like(lma_values[0,:,:])


    # Creating plots for all PFTs for chl
    plt.figure()
    colors = ['g', 'b', 'r']  # Colors for PFTs 2, 3, 4
    legend_labels = []

    for idx, pft_category in enumerate(pft_categories):
        mask = (pft_values == pft_category)[:, :]
        new_chl_values[mask] = average_chl_per_pft[idx]
        new_std_chl_values[mask] = std_chl_per_pft[idx]
        #new_lai_values[mask] = average_lai_per_pft[idx]
        #new_std_lai_values[mask] = std_lai_per_pft[idx]
        #new_lma_values[mask] = average_lma_per_pft[idx]
        #new_std_lma_values[mask] = std_lma_per_pft[idx]

        # Extract chl values for the current PFT category
        chl_for_pft = masked_chl_values[pft_values == pft_category].compressed()
        avg_chl = average_chl_per_pft[idx]
        std_chl = std_chl_per_pft[idx]

        # Create a histogram for chl values with transparency
        plt.hist(chl_for_pft, bins=30, density=True, alpha=0.5, color=colors[idx])

        # Add a vertical line for the average
        plt.axvline(avg_chl, color=colors[idx], linestyle='dashed', linewidth=1)

        # Prepare text for the legend
        avg_std_text = f'{avg_chl:.2f} ± {std_chl:.2f}'
        legend_labels.append(f'PFT {pft_category}: {avg_std_text}')

    plt.title(f'Combined CHL Distribution ({dates_without_times[int(time)]})')
    plt.xlabel(r'CHL Value ($\mu$g.cm$^{-2}$)')
    plt.ylabel('Density')
    plt.ylim(0, 0.5)  # Set y-axis limits
    plt.xlim(0, 20)  # Set y-axis limits
    plt.legend(legend_labels)

    # Save the combined plot as a PNG file
    plt.savefig(f'combined_chl_histogram_time_{time}.png')
    plt.close()

    # Creating plots for all PFTs for lai
    plt.figure()
    colors = ['g', 'b', 'r']  # Colors for PFTs 2, 3, 4
    legend_labels = []

    for idx, pft_category in enumerate(pft_categories):
        mask = (pft_values == pft_category)[:, :]
        #new_chl_values[mask] = average_chl_per_pft[idx]
        #new_std_chl_values[mask] = std_chl_per_pft[idx]
        new_lai_values[mask] = average_lai_per_pft[idx]
        new_std_lai_values[mask] = std_lai_per_pft[idx]
        #new_lma_values[mask] = average_lma_per_pft[idx]
        #new_std_lma_values[mask] = std_lma_per_pft[idx]

        # Extract chl values for the current PFT category
        lai_for_pft = masked_lai_values[pft_values == pft_category].compressed()
        avg_lai = average_lai_per_pft[idx]
        std_lai = std_lai_per_pft[idx]

        # Create a histogram for chl values with transparency
        plt.hist(lai_for_pft, bins=30, density=True, alpha=0.5, color=colors[idx])

        # Add a vertical line for the average
        plt.axvline(avg_lai, color=colors[idx], linestyle='dashed', linewidth=1)

        # Prepare text for the legend
        avg_std_text = f'{avg_lai:.2f} ± {std_lai:.2f}'
        legend_labels.append(f'PFT {pft_category}: {avg_std_text}')

    plt.title(f'Combined LAI Distribution ({dates_without_times[int(time)]})')
    plt.xlabel(r'LAI Value (m$^{2}$.m$^{-2}$)')
    plt.ylabel('Density')
    plt.ylim(0, 1.0)  # Set y-axis limits
    plt.xlim(0, 7.5)  # Set y-axis limits
    plt.legend(legend_labels)

    # Save the combined plot as a PNG file
    plt.savefig(f'combined_lai_histogram_time_{time}.png')
    plt.close()

    # Creating plots for all PFTs for lma
    plt.figure()
    colors = ['g', 'b', 'r']  # Colors for PFTs 2, 3, 4
    legend_labels = []

    for idx, pft_category in enumerate(pft_categories):
        mask = (pft_values == pft_category)[:, :]
        #new_chl_values[mask] = average_chl_per_pft[idx]
        #new_std_chl_values[mask] = std_chl_per_pft[idx]
        #new_lai_values[mask] = average_lai_per_pft[idx]
        #new_std_lai_values[mask] = std_lai_per_pft[idx]
        new_lma_values[mask] = average_lma_per_pft[idx]
        new_std_lma_values[mask] = std_lma_per_pft[idx]

        # Extract chl values for the current PFT category
        lma_for_pft = masked_lma_values[pft_values == pft_category].compressed()
        avg_lma = average_lma_per_pft[idx]
        std_lma = std_lma_per_pft[idx]

        # Create a histogram for chl values with transparency
        plt.hist(lma_for_pft, bins=30, density=True, alpha=0.5, color=colors[idx])

        # Add a vertical line for the average
        plt.axvline(avg_lma, color=colors[idx], linestyle='dashed', linewidth=1)

        # Prepare text for the legend
        avg_std_text = f'{avg_lma:.2e} ± {std_lma:.2e}'
        legend_labels.append(f'PFT {pft_category}: {avg_std_text}')

    print(np.min(lma_for_pft),np.max(lma_for_pft))
    plt.title(f'Combined LMA Distribution ({dates_without_times[int(time)]})')
    plt.xlabel(r'LMA Value (g.cm$^{-2}$)')
    plt.ylabel('Density')
    plt.ylim(0, 800)  # Set y-axis limits
    plt.xlim(0, 0.03)  # Set y-axis limits
    plt.legend(legend_labels)

    # Save the combined plot as a PNG file
    plt.savefig(f'combined_lma_histogram_time_{time}.png')
    plt.close()

    

    # Create new xarray DataArrays for chl and std
    new_chl_da = xr.DataArray(new_chl_values[np.newaxis,:,:], dims=('time', 'lat', 'lon'), coords={'time': chl_ds['time'], 'lat': chl_ds['lat'], 'lon': chl_ds['lon']})
    new_std_chl_da = xr.DataArray(new_std_chl_values[np.newaxis,:,:], dims=('time', 'lat', 'lon'), coords={'time': chl_ds['time'], 'lat': chl_ds['lat'], 'lon': chl_ds['lon']})


    # Replace zeros by NaN in new maps
    new_chl_da = new_chl_da.where(new_chl_da != 0, np.nan)
    new_std_chl_da = new_std_chl_da.where(new_std_chl_da != 0, np.nan)

    # Create a new xarray Dataset with chl and std DataArrays
    new_chl_data = {'chl': new_chl_da, 'std': new_std_chl_da}
    new_chl_ds = xr.Dataset(new_chl_data)
    

    # Save the new chl map with masked values and standard deviation to a NetCDF file
    new_chl_ds.to_netcdf(f'masked_chl_aviris_dangermond_time_{time}.nc')

    print(lai_values.shape)
    print(new_lai_values.shape)
    # Create new xarray DataArrays for lai and std
    #new_lai_da = xr.DataArray(new_lai_values[np.newaxis,:,:], dims=('time', 'lat', 'lon'), coords={'time': lai_ds['time'], 'lat': lai_ds['lat'], 'lon': lai_ds['lon']})
    lai_values = np.clip(lai_values, 1e-1, 20.)  # Clip values between 0.1 and max_lai
    new_lai_values = lai_values[0,:,:]*(new_std_lai_values/new_std_lai_values)
    
    new_lai_da = xr.DataArray(new_lai_values[np.newaxis,:,:], dims=('time', 'lat', 'lon'), coords={'time': lai_ds['time'], 'lat': lai_ds['lat'], 'lon': lai_ds['lon']})
    new_std_lai_da = xr.DataArray(new_std_lai_values[np.newaxis,:,:], dims=('time', 'lat', 'lon'), coords={'time': lai_ds['time'], 'lat': lai_ds['lat'], 'lon': lai_ds['lon']})
    
    # Replace zeros by NaN in new maps
    new_lai_da = new_lai_da.where(new_lai_da != 0, np.nan)
    new_lai_da = new_lai_da.where(~np.isnan(new_lai_da), np.nan)

    new_std_lai_da = new_std_lai_da.where(new_std_lai_da != 0, np.nan)

    # Create a new xarray Dataset with lai and std DataArrays
    new_lai_data = {'lai': new_lai_da, 'std': new_std_lai_da}
    new_lai_ds = xr.Dataset(new_lai_data)
    
    # Save the new lai map with masked values and standard deviation to a NetCDF file
    new_lai_ds.to_netcdf(f'masked_lai_aviris_dangermond_time_{time}_v1.nc')

    # Create new xarray DataArrays for lma and std
    new_lma_da = xr.DataArray(new_lma_values[np.newaxis,:,:], dims=('time', 'lat', 'lon'), coords={'time': lma_ds['time'], 'lat': lma_ds['lat'], 'lon': lma_ds['lon']})
    new_std_lma_da = xr.DataArray(new_std_lma_values[np.newaxis,:,:], dims=('time', 'lat', 'lon'), coords={'time': lma_ds['time'], 'lat': lma_ds['lat'], 'lon': lma_ds['lon']})
    
    new_lma_da = new_lma_da.where(new_lma_da != 0, np.nan)
    new_std_lma_da = new_std_lma_da.where(new_std_lma_da != 0, np.nan)

    # Create a new xarray Dataset with lma and std DataArrays
    new_lma_data = {'lma': new_lma_da, 'std': new_std_lma_da}
    new_lma_ds = xr.Dataset(new_lma_data)
    


    # Save the new lma map with masked values and standard deviation to a NetCDF file
    new_lma_ds.to_netcdf(f'masked_lma_aviris_dangermond_time_{time}.nc')
    
    
    # Informative message: Data processing and saving completed for the current time
    print(f'Data processing and saving completed for time: {time}\n')

    # Close the opened datasets
    chl_ds.close()
    lai_ds.close()
    lma_ds.close()
    pft_ds.close()

