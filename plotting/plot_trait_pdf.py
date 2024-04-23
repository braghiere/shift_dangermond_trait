import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Define input file paths
pft_file = '/net/fluo/data1/students/renato/aviris_dangermond/California_Vegetation_WHRTYPE_Dangermond/output_latlon.nc'
chlorophyll_file = 'chl_aviris_dangermond_time_00.nc'
output_pdf = 'chlorophyll_per_pft.pdf'

# Open the PFT file using xarray
pft_ds = xr.open_dataset(pft_file)

# Open the chlorophyll file using xarray
chl_ds = xr.open_dataset(chlorophyll_file)

# Extract latitude and longitude values
lat_values = chl_ds['latitude'].values
lon_values = chl_ds['longitude'].values

# Initialize a PDF for plots
pdf_pages = PdfPages(output_pdf)

# Iterate over PFT values and create plots
for pft_index in range(pft_ds['Band1'].shape[0]):
    pft_name = pft_ds['Band1'][pft_index].values.astype(str)
    chl_values = chl_ds['chl'].values[:, :, :]
    chl_mean = np.nanmean(chl_values, axis=(1, 2))
    chl_std = np.nanstd(chl_values, axis=(1, 2))

    # Create a plot for chlorophyll values per PFT
    plt.figure(figsize=(10, 6))
    plt.plot(chl_mean, label='Mean Chlorophyll')
    plt.fill_between(range(len(chl_mean)), chl_mean - chl_std, chl_mean + chl_std, alpha=0.3, label='Std Dev')
    plt.title(f'Chlorophyll per PFT: {pft_name}')
    plt.xlabel('Time')
    plt.ylabel('Chlorophyll (micromol.m-2)')
    plt.legend()
    
    # Save the plot to the PDF
    pdf_pages.savefig()
    plt.close()

# Close the PDF file
pdf_pages.close()

print(f'Chlorophyll per PFT PDF saved to {output_pdf}')

