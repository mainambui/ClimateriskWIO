#!/bin/bash

# Define the input directory and start and end years
input_dirs=("CDD" "EVAP" "NPP" "Rainfall" "SST" "Temperature" "tx90")
start_year="2021"
end_year="2050"

# Loop through the input directories
for input_dir in "${input_dirs[@]}"; do
    # Loop through the files in the input directory with the extension .nc
    for file in "${input_dir}"/*.nc; do
        # Extract the file name without the path and extension
        file_name=$(basename "$file" .nc)
        # Define the output file name
        output_file="1_ESMs/${start_year}-${end_year}/${file_name}.nc"
        # Run the cdo command with the input and output file names
        cdo -yearmean "$file" "$output_file"
    done
done
