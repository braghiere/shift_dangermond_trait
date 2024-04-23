import xarray as xr
import pandas as pd
from sklearn.metrics import mean_squared_error, r2_score
import numpy as np

# Open the NetCDF files using xarray
file1 = "/net/fluo/data3/data/FluoData1/students/renato/aviris_dangermond/compare_sif/shift_fluxes_day_00_reg_tropomi.nc"
file2 = "/net/fluo/data3/data/FluoData1/students/renato/aviris_dangermond/TROPOMI_dangermond/TROPOMI_SIF740nm-v1.001deg_regrid_Dangermond_tll_clipped.nc"

# Open datasets
ds1 = xr.open_dataset(file1, decode_times=False)  # Use decode_times=False to handle time variable as numeric
ds2 = xr.open_dataset(file2, decode_times=False)

# Convert target date to numeric format compatible with the 'time' variable in ds2
target_date = pd.Timestamp("2022-02-24").value / 10**9  # Convert to seconds since epoch

# Select the data from the dataset with multiple times based on the specific date
ds2_selected = ds2.sel(time=target_date, method="nearest")
ds1_selected = ds1['sif740'].squeeze(dim='time')

# Filter out NaN values
valid_indices = ~np.isnan(ds1_selected.values) & ~np.isnan(ds2_selected['sif'].values)
#observed = ds2_selected['sif'].values[valid_indices]
observed = ds2_selected['sif_relative'].values[valid_indices]
predicted = ds1_selected.values[valid_indices]

# Check if there are still valid values after filtering
if len(observed) == 0 or len(predicted) == 0:
    print("No valid data points for comparison.")
else:
    # Calculate mean absolute difference
    mean_abs_diff = np.mean(np.abs(predicted - observed))
    print("Mean Absolute Difference:", mean_abs_diff)

    # Calculate R^2 score
    r2 = r2_score(observed, predicted)
    print("R^2 Score:", r2)



import matplotlib.pyplot as plt

# Create a scatter plot
plt.figure(figsize=(8, 6))
plt.scatter(observed, predicted, color='blue', label='Observed vs Predicted')
plt.plot(observed, observed, color='red', linestyle='--', label='Perfect Fit')

plt.title('Observed vs Predicted SIF Values')
plt.xlabel('Observed SIF')
plt.ylabel('Predicted SIF')
plt.legend()
plt.grid(True)
plt.show()
