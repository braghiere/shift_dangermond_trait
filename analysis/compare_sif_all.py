import xarray as xr
import pandas as pd
from sklearn.metrics import mean_squared_error, r2_score
import numpy as np
import matplotlib.pyplot as plt

# Given dates
dates = ["2022-02-24T00:00:00.000000", "2022-02-28T00:00:00.000000", "2022-03-08T00:00:00.000000",
         "2022-03-16T00:00:00.000000", "2022-03-22T00:00:00.000000", "2022-04-05T00:00:00.000000",
         "2022-04-12T00:00:00.000000", "2022-04-20T00:00:00.000000", "2022-04-29T00:00:00.000000",
         "2022-05-03T00:00:00.000000", "2022-05-11T00:00:00.000000", "2022-05-17T00:00:00.000000",
         "2022-05-29T00:00:00.000000"]

# Base file path
base_file_path = "/net/fluo/data3/data/FluoData1/students/renato/aviris_dangermond/compare_sif/shift_fluxes_day_{}_reg_tropomi.nc"

# Open the NetCDF file for TROPOMI dataset
file2 = "/net/fluo/data3/data/FluoData1/students/renato/aviris_dangermond/TROPOMI_dangermond/TROPOMI_SIF740nm-v1.001deg_regrid_Dangermond_tll_clipped.nc"
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

# Calculate mean absolute difference and R^2 score for all accumulated data points
mean_abs_diff = np.mean(np.abs(np.array(predicted_all) - np.array(observed_all)))
r2 = r2_score(np.array(observed_all), np.array(predicted_all))

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
plt.show()

