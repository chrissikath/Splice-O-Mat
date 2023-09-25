#!/bin/bash

# Find all folders in the "star/" directory
for folder in star/*/; do
    # Check if the file "Aligned.out.bam" exists in the folder
    if [ -f "$folder/Aligned.out.bam" ]; then
        # Delete the file "Aligned.out.bam"
        rm "$folder/Aligned.out.bam"
        echo "Deleted Aligned.out.bam in $folder"
    else
        echo "Aligned.out.bam not found in $folder"
    fi
done

