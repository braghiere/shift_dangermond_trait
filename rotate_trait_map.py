import xarray as xr

# Function to transpose and save the dataset
def transpose_netcdf(file_path, output_path):
    ds = xr.open_dataset(file_path)
    transposed_ds = ds.transpose('lat', 'lon')
    transposed_ds.to_netcdf(output_path)

# Loop through times 00 to 11
for time in range(0,13):
    # Properly format the time with two digits
    time_str = str(time).zfill(2)
    input_file = f'shift_fittings_lwc_time_{time_str}.nc'
    output_file = f'transposed_shift_fittings_lwc_time_{time_str}.nc'
    transpose_netcdf(input_file, output_file)
    print(f'Processed {input_file} -> {output_file}')

