#!/bin/bash

# Enhance lat/lon definitions in the original NetCDF file
#for file in shift_fluxes_day_*_reg_jmax.nc; do
for file in pft_shift_fluxes_day_*_jmax.nc; do
    ncap2 -s 'lat@units="degrees_north";lat@standard_name="latitude";lon@units="degrees_east";lon@standard_name="longitude";' $file ${file%.nc}_fixed.nc
done

# Define the target grid file from TROPOMI data
TARGET_GRID_FILE="target_grid_TROPOMI.txt"
TARGET_NETCDF_FILE="/net/fluo/data3/data/FluoData1/students/renato/aviris_dangermond/TROPOMI_dangermond/TROPOMI_SIF740nm-v1.001deg_regrid_Dangermond_tll_clipped.nc"

# Generate a grid description file from the target NetCDF
cdo griddes "$TARGET_NETCDF_FILE" > "$TARGET_GRID_FILE"

# Loop through each 'shift_fluxes' and 'pft_shift' file and regrid it
#for file in shift_fluxes_day_*_reg_jmax_fixed.nc; do
for file in pft_shift_fluxes_day_*_jmax_fixed.nc; do
    if [ -f "$file" ]; then
        # Perform the remapping
        OUTPUT_FILE="${file%.nc}_regridded.nc"
        cdo remapbil,"$TARGET_GRID_FILE" "$file" "$OUTPUT_FILE"
        echo "Regridded $file to $OUTPUT_FILE"
    else
        echo "File $file does not exist"
    fi
done


