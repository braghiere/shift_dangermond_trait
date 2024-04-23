#!/bin/bash

# Create folders
for i in {0..7}; do
    folder_name=$(printf "%02d_%02d" $((i*10)) $((i*10+9)))
    mkdir "$folder_name"
done

# Move files to respective folders
for file in *time_01_??_??.nc; do
    if [[ -f $file ]]; then
        num=$(echo "$file" | grep -oP 'time_01_\K\d{2}(?=_\d{2}_\d{2}\.nc)')
        folder_name=$(printf "%02d_%02d" $((num/10*10)) $(((num/10*10)+9)))
        mv "$file" "$folder_name/"
    fi
done

echo "Files moved to respective folders."

