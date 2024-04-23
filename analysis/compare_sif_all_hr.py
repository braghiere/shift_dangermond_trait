import xarray as xr
import pandas as pd
from sklearn.metrics import mean_squared_error, r2_score
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.stats import pearsonr

import sys

# Given dates
dates = ["2022-02-24T00:00:00.000000", "2022-02-28T00:00:00.000000", "2022-03-08T00:00:00.000000",
         "2022-03-16T00:00:00.000000", "2022-03-22T00:00:00.000000", "2022-04-05T00:00:00.000000",
         "2022-04-12T00:00:00.000000", "2022-04-20T00:00:00.000000", "2022-04-29T00:00:00.000000",
         "2022-05-03T00:00:00.000000", "2022-05-11T00:00:00.000000", "2022-05-17T00:00:00.000000",
         "2022-05-29T00:00:00.000000"]

# Base file path
base_file_path = "/net/fluo/data3/data/FluoData1/students/renato/aviris_dangermond/fitting/shift_fluxes_day_{}_reg.nc"

# Open the NetCDF file for TROPOMI dataset
file2 = "/net/fluo/data3/data/FluoData1/students/renato/aviris_dangermond/TROPOMI_dangermond/TROPOMI_SIF740nm-v1.001deg_regrid_Dangermond_tll_clipped_458_492.nc"
ds2 = xr.open_dataset(file2, decode_times=False)

# Initialize empty lists to accumulate data points
observed_all = []
predicted_all = []

for i, target_date_str in enumerate(dates):
    # Construct file1 dynamically based on the time step index (i)
    file1 = base_file_path.format(str(i).zfill(2))

    # Open dataset for file1
    ds1 = xr.open_dataset(file1, decode_times=False)  
    
    target_date = pd.Timestamp(target_date_str).value / 10**9  # Convert to seconds since epoch

    # Select the data from the dataset with multiple times based on the specific date
    ds2_selected = ds2.sel(time=target_date, method="nearest")
    ds1_selected = ds1['sif740']
    
     # Filter out NaN values
    valid_indices = ~np.isnan(ds2_selected['sif'].values.flatten()) & ~np.isnan(ds1_selected.values.flatten())
    observed = ds2_selected['sif'].values.flatten()[valid_indices]
    predicted = ds1_selected.values.flatten()[valid_indices]

    # Calculate spatial differences
    spatial_diff_num = observed - predicted

    # Squeeze the array to remove the singleton dimension
    spatial_diff_num = np.squeeze(spatial_diff_num)   
    
    # Calculate spatial differences
    spatial_diff = ds2_selected['sif'].values - ds1_selected.values
    
    # Squeeze the array to remove the singleton dimension
    spatial_diff = np.squeeze(spatial_diff)

    #Get latitude and longitude values
    lat_values = ds2_selected['lat'].values
    lon_values = ds2_selected['lon'].values
    
    # Calculate metrics
    bias = np.mean(spatial_diff_num)
    rmse = np.sqrt(mean_squared_error(observed, predicted))
    #r2 = r2_score(observed, predicted)
    # Calculate Pearson correlation coefficient
    pearson_corr, _ = pearsonr(observed, predicted)

    # Calculate R^2
    r2 = pearson_corr**2



    # Plot spatial differences with latitude and longitude values on the y and x axes
    plt.figure(figsize=(8, 6))
    plt.imshow(spatial_diff, cmap='coolwarm', vmin=-2, vmax=2)  # Adjust vmin and vmax according to your data range
    plt.colorbar()
    plt.title(f'Spatial Differences at {target_date_str}')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    
    
    # Set 2 decimal places for latitude and longitude ticks
    lat_formatter = ticker.FuncFormatter(lambda x, pos: '{:.2f}'.format(lat_values[x]))
    lon_formatter = ticker.FuncFormatter(lambda x, pos: '{:.2f}'.format(lon_values[x]))

    
    # Set 4 ticks for latitude and longitude
    num_ticks = 4
    lat_indices = np.linspace(0, len(lat_values) - 1, num_ticks, dtype=int)
    lon_indices = np.linspace(0, len(lon_values) - 1, num_ticks, dtype=int)

    # Set x and y ticks with actual lat and lon values formatted to 2 decimal places
    plt.xticks(lon_indices, lon_values[lon_indices])
    plt.yticks(lat_indices, lat_values[lat_indices])

    # Set formatted latitude and longitude tick labels
    plt.gca().get_xaxis().set_major_formatter(lon_formatter)
    plt.gca().get_yaxis().set_major_formatter(lat_formatter)

    # Add text box with statistics in the upper right corner
    bbox_props = dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="black", linewidth=0.5)
    plt.gca().text(0.95, 0.95, f'Bias: {bias:.2f}\nRMSE: {rmse:.2f}\nR²: {r2:.2f}',
               transform=plt.gca().transAxes, bbox=bbox_props, verticalalignment='top', horizontalalignment='right')

    
    
    # Save the figure
    plt.savefig(f'spatial_differences_{i}.png')
    #plt.show()
    #sys.exit()
    
    # Close the figure to avoid displaying multiple plots at once
    plt.close()
    print('Figure spatial_differences_{i}.png saved!')
    

    # Filter out NaN values
    valid_indices = ~np.isnan(ds1_selected.values) & ~np.isnan(ds2_selected['sif'].values)
    observed = ds2_selected['sif'].values.reshape(-1)[valid_indices.reshape(-1)]
    predicted = ds1_selected.values.reshape(-1)[valid_indices.reshape(-1)]

    # Check if there are still valid values after filtering
    if len(observed) > 0 and len(predicted) > 0:
        # Accumulate observed and predicted values
        observed_all.extend(observed)
        predicted_all.extend(predicted)

    # Close the datasets to free up resources
    ds1.close()
    ds2_selected.close()
    
    print('Total matrix ammended!')

# Calculate mean absolute difference and R^2 score for all accumulated data points
mean_abs_diff = np.mean(np.abs(np.array(predicted_all) - np.array(observed_all)))
#r2 = r2_score(np.array(observed_all), np.array(predicted_all))
# Calculate Pearson correlation coefficient
pearson_corr, _ = pearsonr(np.array(observed_all), np.array(predicted_all))

# Calculate R^2
r2 = pearson_corr**2

print("Mean Absolute Difference:", mean_abs_diff)
print("R^2 Score:", r2)

# Create a scatter plot for all accumulated data points
plt.figure(figsize=(8, 6))
plt.scatter(observed_all, predicted_all, color='blue', label='Observed vs Predicted')
plt.plot(observed_all, observed_all, color='red', linestyle='--', label='Perfect Fit')

plt.title('Observed vs Predicted SIF Values')
plt.xlabel('Observed SIF')
plt.ylabel('Predicted SIF')
plt.legend()
plt.grid(True)
plt.savefig(f'scatter_plot_all_trait.png')
plt.close()

