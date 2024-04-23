#!/bin/bash

# Define the variable names
variables=('lai' 'lma' 'chl')

# Loop through dates and variables
for date_index in {00..12}; do
    for variable in "${variables[@]}"; do
        # Input and output file paths
        input_file="${variable}_aviris_dangermond_time_${date_index}.nc"
        output_file="${variable}_aviris_dangermond_time_${date_index}_reg.nc"
        
        # Run CDO remapbil command
        cdo remapbil,../California_Vegetation_WHRTYPE_Dangermond/output_latlon.nc "$input_file" "$output_file"
        
        echo "Remapped $variable data for 20$(printf "%02d" $date_index)"
    done
done

echo "Remapping process completed!"

