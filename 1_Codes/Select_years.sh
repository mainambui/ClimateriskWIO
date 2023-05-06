#!/bin/bash

# Define the start and end years
start_year=2021
end_year=2050

# Loop through the files
for file in Annual_means/*.nc; do
    # Extract the file name without the path and extension
    file_name=$(basename "$file" .nc)
    # Define the output file name
    output_file="${start_year}-${end_year}/${file_name}.nc"
    # Run the cdo command with the input and output file names
    cdo timmean -selyear,"$start_year"/"$end_year" "$file" "$output_file"
done
