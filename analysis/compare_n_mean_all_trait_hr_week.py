import xarray as xr
import pandas as pd
from sklearn.metrics import mean_squared_error, r2_score
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Polygon
from scipy.stats import pearsonr
import geopandas as gpd
import sys
from mpl_toolkits.basemap import Basemap
import pyproj
from pyproj import Proj, transform
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FuncFormatter

from scipy.stats import linregress
from sklearn.metrics import mean_absolute_error
import statsmodels.api as sm

import scipy.stats as st 


# Define the source and target CRS
src_crs = Proj('EPSG:4326')
target_crs = Proj(proj='latlong', datum='WGS84')

# Now, use lon_values_reprojected and lat_values_reprojected in your plt.imshow() and shapefile plotting functions

# Load the shapefile using geopandas
shapefile_path = "/net/fluo/data3/data/FluoData1/students/renato/aviris_dangermond/dangermond_bound/dan_bound2.shp"
gdf = gpd.read_file(shapefile_path)
#gdf = gdf.to_crs('EPSG:2229')
#gdf = gdf.to_crs('EPSG:4326')

#shape_values = gdf.values

#print(shape_values)

#sys.exit()

#stateshp = gpd.read_file('/net/fluo/data3/data/FluoData1/students/renato/aviris_dangermond/dangermond_bound/dan_bound.shp')
#print stateshp1.crs
# Reproject to EPSG:4326
#stateshp = stateshp.to_crs(epsg=4326)
#stateshp.to_file('/net/fluo/data3/data/FluoData1/students/renato/aviris_dangermond/dangermond_bound/dan_bound1.shp', driver='ESRI Shapefile')

#sys.exit()
# Print the information about the shapefile
#print(gdf.info())

# Given dates
dates = ["2022-02-24T00:00:00.000000", "2022-02-28T00:00:00.000000", "2022-03-08T00:00:00.000000",
         "2022-03-16T00:00:00.000000", "2022-03-22T00:00:00.000000", "2022-04-05T00:00:00.000000",
         "2022-04-12T00:00:00.000000", "2022-04-20T00:00:00.000000", "2022-04-29T00:00:00.000000",
         "2022-05-03T00:00:00.000000", "2022-05-11T00:00:00.000000", "2022-05-17T00:00:00.000000",
         "2022-05-29T00:00:00.000000"]

# Base file path
base_file_path = "/net/fluo/data3/data/FluoData1/students/renato/aviris_dangermond/fitting/shift_fluxes_day_{}_reg_h_lma.nc"

# Open the NetCDF file for TROPOMI dataset
file2 = "/net/fluo/data3/data/FluoData1/students/renato/aviris_dangermond/TROPOMI_dangermond/TROPOMI_SIF740nm-v1.001deg_regrid_Dangermond_tll_clipped_458_492.nc"
ds2 = xr.open_dataset(file2, decode_times=True)

# Save dates, mean, and std values to a text file
try:
    file = open('mean_n_trait_statistics_hr.txt', 'w')
    file.write('Date\tMean Observed\tStd Observed\tCI95 Observed\tMean Predicted\tStd Predicted\tCI95 Predicted\n')
except Exception as e:
    print(f"Error: {e}")


# Initialize empty lists to accumulate data points
observed_all = []
predicted_all = []

# Initialize a list to store individual spatial differences
spatial_diff_list = []

for i, target_date_str in enumerate(dates):
    # Construct file1 dynamically based on the time step index (i)
    file1 = base_file_path.format(str(i).zfill(2))

    # Open dataset for file1
    ds1 = xr.open_dataset(file1, decode_times=False)  
    
    target_date = pd.Timestamp(target_date_str).value / 10**9  # Convert to seconds since epoch
    print('target_date=',target_date)
    # Convert target_date to datetime64[ns]
    target_date = pd.to_datetime(target_date, unit='s')
    print('target_date=',target_date)

    # Define a time range for a week centered around the target date
    start_date = target_date - pd.Timedelta(days=3)  # 3 days before the target date
    end_date = target_date + pd.Timedelta(days=3)    # 3 days after the target date

    # Select the data from the dataset with multiple times based on the specific date
    ds2_selected = ds2.sel(time=target_date, method="nearest")
    print(ds2_selected)
    
    # Select the data from ds2 within the specified time range
    ds2_selected = ds2.sel(time=slice(start_date, end_date))
    ds2_selected = ds2_selected.mean(dim='time', skipna=True)
    print(ds2_selected)
    #sys.exit()

    # Calculate the average time of the selected data in ds2_selected

    ds1_selected = ds1['sif740']
    
    # Filter out NaN values
    valid_indices = ~np.isnan(ds2_selected['n'].values.flatten()) & ~np.isnan(ds1_selected.values.flatten())
    if valid_indices.any():
        observed = ds2_selected['n'].values.flatten()[valid_indices]
        predicted = ds1_selected.values.flatten()[valid_indices]
        
        print("Observed shape before cleaning:", observed.shape)
        print("Predicted shape before cleaning:", predicted.shape)

        observed = observed[~np.isnan(observed)]
        predicted = predicted[~np.isnan(predicted)]

        print("Observed shape after cleaning:", observed.shape)
        print("Predicted shape after cleaning:", predicted.shape)


        print('mean obs =', np.mean(observed))
        print('mean pred =', np.mean(predicted))

        # Calculate spatial differences
        spatial_diff_num =  observed 

        # Squeeze the array to remove the singleton dimension
        spatial_diff_num = np.squeeze(spatial_diff_num)   
        
        # Calculate spatial differences
        spatial_diff =  ds2_selected['n'].values
        
        # Squeeze the array to remove the singleton dimension
        spatial_diff = np.squeeze(spatial_diff)
        
        
        # Calculate spatial differences and append them to the list
        spatial_diff_list.append(spatial_diff)

        #Get latitude and longitude values
        lat_values = ds2_selected['lat'].values
        lon_values = ds2_selected['lon'].values
        
        # Calculate metrics
        mae = mean_absolute_error(observed, predicted)

        # Calculate mean, std, bias, RMSE, and R²
        mean_observed = np.mean(observed)
        mean_predicted = np.mean(predicted)
        std_observed = np.std(observed)
        std_predicted = np.std(predicted)
        bias = mean_predicted - mean_observed
        rmse = np.sqrt(np.mean((predicted - observed)**2))
        slope, intercept, r_value, p_value, std_err = linregress(observed, predicted)
        r2 = r_value**2


        
        #plt.imshow(spatial_diff)
        #plt.colorbar()
        #plt.show()
        #plt.close()
        print("mae:", mae)
        print("bias:", bias)
        print("RMSE:", rmse)
        print("R^2 Score:", r2)

        # create 95% confidence interval 
        ci95_observed = st.t.interval(0.95, df=len(observed)-1, 
                        loc=np.mean(observed), 
                        scale=st.sem(observed)) 

        # create 95% confidence interval 
        ci95_predicted = st.t.interval(0.95, df=len(predicted)-1, 
                        loc=np.mean(predicted), 
                        scale=st.sem(predicted)) 
  
        # Append statistics for the current date to the file
        file.write(f'{target_date_str}\t{mean_observed}\t{std_observed}\t{ci95_observed}\t{mean_predicted}\t{std_predicted}\t{ci95_predicted}\t{bias}\t{rmse}\t{r2}\n')

        # Create a Basemap instance with the desired projection
        m = Basemap(projection='merc', llcrnrlat=lat_values.min(), urcrnrlat=lat_values.max(),
                llcrnrlon=lon_values.min(), urcrnrlon=lon_values.max(), resolution='c')

        x,y = np.meshgrid(lon_values, lat_values) 
        X,Y = m(x, y)
        
        # Define the number of decimal places for meridian and parallel labels
        decimal_places = 3

        # Define a custom formatter function to format the labels with the specified number of decimal places
        def format_labels(x, pos):
            return '{:.{decimal_places}f}'.format(x, decimal_places=decimal_places)


        # Create a figure and axes
        plt.figure(figsize=(8, 6))
        m.drawparallels(np.arange(-90.,91.,0.025), labels=[1,0,0,1],    dashes=[1,1], linewidth=0.25, color='0.5',fontsize=8)
        m.drawmeridians(np.arange(0., 360., 0.025), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5',fontsize=8)
        #m.drawcoastlines(color='0.6', linewidth=0.5)

        # Plot spatial differences using pcolormesh with Basemap
        cs = m.pcolormesh(X, Y, spatial_diff, cmap='viridis', vmin=0, vmax=0.05)

        vmin = 0.
        vmax = 0.05
        cbar = m.colorbar(cs, pad='10%',ticks=np.linspace(vmin,vmax,7),format='%.2f')
        cbar.ax.get_yaxis().labelpad = 10
        cbar.ax.set_ylabel(r'N', rotation=270, verticalalignment='center', color='black', size=16)
        cbar.solids.set_edgecolor("face")
        #cbar.set_clim(vmin,vmax)
        # Assuming cbar is your Colorbar object
        #cbar.ax.set_clim(vmin, vmax)
        cbar.ax.tick_params(labelsize='large')
        #gdf.plot(ax=m, linewidth=1, edgecolor='black', facecolor='none', legend=True)
        #m.readshapefile(shapefile_path,'Geometry')
        # Plot contour of Dangermond shapefile
        #m.readshapefile(shapefile_path, 'dan_bound2', linewidth=2, color='black')
        
        # Manually draw shapefile polygons on top of the contour plot
        for shape in gdf['geometry']:
            if shape.geom_type == 'Polygon':
                x, y = m(shape.exterior.coords.xy[0], shape.exterior.coords.xy[1])
                polygon = Polygon(list(zip(x, y)), edgecolor='black', linewidth=1.5, facecolor='none')
                plt.gca().add_patch(polygon)

        # Plot spatial differences with latitude and longitude values on the y and x axes
        #plt.imshow(spatial_diff, cmap='coolwarm', vmin=-1, vmax=1)  # Adjust vmin and vmax according to your data range
    
        # Overlay shapefile contour on top of the spatial differences plot
        #gdf.plot(ax=plt.gca(), linewidth=1, edgecolor='black', facecolor='none', legend=True)


        # Set title, labels, and tick formatters
        plt.title(f'N at {target_date_str}',fontsize=14)
        #plt.xlabel('Longitude')
        #plt.ylabel('Latitude')
    

        # Add colorbar to the right of the plot
        #divider = make_axes_locatable(ax)
        #cax = divider.append_axes("right", size="5%", pad=0.1)


        
        # Set 2 decimal places for latitude and longitude ticks
        #lat_formatter = ticker.FuncFormatter(lambda x, pos: '{:.2f}'.format(lat_values[x]))
        #lon_formatter = ticker.FuncFormatter(lambda x, pos: '{:.2f}'.format(lon_values[x]))
        
        
        # Get latitude and longitude values from your SIF data
        #lat_min, lat_max = lat_values.min(), lat_values.max()
        #lon_min, lon_max = lon_values.min(), lon_values.max()
        #print( lon_min, lon_max,lat_min, lat_max)
        
        # Set 4 ticks for latitude and longitude
        #num_ticks = 4
        #lat_indices = np.linspace(0, len(lat_values) - 1, num_ticks, dtype=int)
        #lon_indices = np.linspace(0, len(lon_values) - 1, num_ticks, dtype=int)

        # Set x and y ticks with actual lat and lon values formatted to 2 decimal places
        #plt.xticks(lon_indices, lon_values[lon_indices])
        #plt.yticks(lat_indices, lat_values[lat_indices])

        # Set formatted latitude and longitude tick labels
        #plt.gca().get_xaxis().set_major_formatter(lon_formatter)
        #plt.gca().get_yaxis().set_major_formatter(lat_formatter)
        


        # Add text box with statistics in the upper right corner
        bbox_props = dict(boxstyle="round,pad=0.8", facecolor="white", edgecolor="black", linewidth=0.5)
        plt.gca().text(0.95, 0.95, f'Mean: {mean_observed:.2f}\nSTD: {std_observed:.2f}',
                transform=plt.gca().transAxes, bbox=bbox_props, verticalalignment='top', horizontalalignment='right',fontsize=14)
                


        # Invert the latitude axis
        #plt.gca().invert_yaxis()
        plt.tight_layout()
        # Save the figure
        plt.savefig(f'n_{i}_hr.png',dpi=300)
        #plt.show()
        #sys.exit()
        
        # Close the figure to avoid displaying multiple plots at once
        plt.close()
        print('Figure spatial_differences_{i}.png saved!')
        

        # Filter out NaN values
        valid_indices = ~np.isnan(ds1_selected.values) & ~np.isnan(ds2_selected['n'].values)
        observed = ds2_selected['n'].values.reshape(-1)[valid_indices.reshape(-1)]
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
        print(spatial_diff_list)
    
    
# Calculate the average spatial difference after the loop
average_spatial_diff = np.nanmean(spatial_diff_list, axis=0)  

# Create a figure and axes
plt.figure(figsize=(8, 6))
m.drawparallels(np.arange(-90., 91., 0.025), labels=[1, 0, 0, 1], dashes=[1, 1], linewidth=0.25, color='0.5', fontsize=8)
m.drawmeridians(np.arange(0., 360., 0.025), labels=[1, 0, 0, 1], dashes=[1, 1], linewidth=0.25, color='0.5', fontsize=8)

# Plot spatial differences using pcolormesh with Basemap
cs = m.pcolormesh(X, Y, average_spatial_diff, cmap='viridis', vmin=0., vmax=0.05)

vmin = 0.
vmax = 0.05
cbar = m.colorbar(cs, pad='10%', ticks=np.linspace(vmin, vmax, 7), format='%.2f')
cbar.ax.get_yaxis().labelpad = 10
cbar.ax.set_ylabel(r'N', rotation=270, verticalalignment='center', color='black', size=16)
cbar.solids.set_edgecolor("face")
cbar.ax.tick_params(labelsize='large')

# Manually draw shapefile polygons on top of the contour plot
for shape in gdf['geometry']:
    if shape.geom_type == 'Polygon':
       x, y = m(shape.exterior.coords.xy[0], shape.exterior.coords.xy[1])
       polygon = Polygon(list(zip(x, y)), edgecolor='black', linewidth=1.5, facecolor='none')
       plt.gca().add_patch(polygon)

# Set title, labels, and tick formatters
plt.title(f'Mean Spatial Difference', fontsize=14)

# Add text box with statistics in the upper right corner
bbox_props = dict(boxstyle="round,pad=0.8", facecolor="white", edgecolor="black", linewidth=0.5)
plt.gca().text(0.95, 0.95, f'Mean: {mean_observed:.2f}\nSTD: {std_observed:.2f}',
           transform=plt.gca().transAxes, bbox=bbox_props, verticalalignment='top', horizontalalignment='right', fontsize=14)


plt.tight_layout()
# Save the figure
plt.savefig(f'n_mean_hr.png', dpi=300)


# Close the figure to avoid displaying multiple plots at once
plt.close()
print(f'Figure pft_spatial_differences_{i}.png saved!')
  
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

# Perform linear regression
slope, intercept, r_value, p_value, std_err = linregress(observed_all, predicted_all)

# Assuming observed_all is a list
observed_all = np.array(observed_all)
predicted_all = np.array(predicted_all)

# Calculate fitted values using the regression equation
fitted_values = intercept + slope * observed_all

# Calculate MAE, bias, RMSE, and R²
mae = mean_absolute_error(observed_all, predicted_all)
bias = np.mean(predicted_all - observed_all)
rmse = np.sqrt(np.mean((predicted_all - observed_all)**2))
r2 = r_value**2

# Fit the linear regression model using statsmodels for AIC and BIC
X = sm.add_constant(observed_all)  # Add a constant term to the predictor
model = sm.OLS(predicted_all, X).fit()

# Get AIC and BIC
aic = model.aic
bic = model.bic

print(aic,bic)

# Plot the regression line
plt.plot(observed_all, fitted_values, color='green', label='Fitted Line')

# Display bias, RMSE, and R² on the plot
plt.text(0.1, 0.9, f'Bias: {bias:.2f}\nRMSE: {rmse:.2f}\nR²: {r2:.2f}\nMAE: {mae:.2f}\nAIC: {aic:.2f}\nBIC: {bic:.2f}', transform=plt.gca().transAxes,
         bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.8'))

print(f'Bias: {bias:.2f}\nRMSE: {rmse:.2f}\nR²: {r2:.2f}\nMAE: {mae:.2f}\nAIC: {aic:.2f}\nBIC: {bic:.2f}')

plt.title('Observed vs Predicted SIF Values')
plt.xlabel('Observed SIF')
plt.ylabel('Predicted SIF')
plt.legend()
plt.grid(True)
#plt.savefig(f'scatter_plot_all_trait_hr.png')
plt.close()


file.close()


