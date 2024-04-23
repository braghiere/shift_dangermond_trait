#!/bin/bash

# Activate Conda environment with cdo
source activate cdo

# Apply cdo command within each folder
for folder in ??_??; do
    if [[ -d "$folder" ]]; then
        cd "$folder"
        folder_number=$(echo "$folder" | sed 's/_//g')

        output_file="${folder_number}_output.nc"
        cdo collgrid,10 *.nc ../"$output_file"
        echo "Command executed in $folder: cdo collgrid,10 $input_file $output_file"

        cd ..
    fi
done

cdo collgrid,8 *.nc output_merged.nc

echo "Files merged!"

