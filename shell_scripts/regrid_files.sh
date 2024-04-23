#!/bin/bash

# Loop through files from 00 to 12
for i in {00..12}; do
    input_file="/net/fluo/data3/data/FluoData1/students/renato/aviris_dangermond/fitting/pft_shift_fluxes_day_${i}_h_lma.nc"
    #input_file="/net/fluo/data3/data/FluoData1/students/renato/aviris_dangermond/fitting/shift_fluxes_day_${i}_reg_h_lma.nc"
    output_file="/net/fluo/data3/data/FluoData1/students/renato/aviris_dangermond/compare_sif/pft_shift_fluxes_day_${i}_reg_tropomi.nc"
    
    # Convert NaN _FillValue to a normal number
    ncatted -a _FillValue,lat,m,f,1.0e36 ${input_file} ${input_file}_temp.nc
    mv ${input_file}_temp.nc ${input_file}
    
    # Run the regrid command
    echo "Regridding ${input_file}..."
    ncremap -i ${input_file} -d /net/fluo/data3/data/FluoData1/students/renato/aviris_dangermond/TROPOMI_dangermond/TROPOMI_SIF740nm-v1.001deg_regrid_Dangermond_tll_clipped.nc -o ${output_file}
    
    # Check if the regrid command was successful
    if [ $? -eq 0 ]; then
        echo "Regridding successful."
    else
        echo "Regridding failed for ${input_file}."
    fi
done

