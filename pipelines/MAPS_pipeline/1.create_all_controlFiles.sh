#!/usr/bin/bash

lwB=$1
upB=$2

# Loop through the prefixes
for prefix in CLK AKA AKY PRD CRP MEL SAN POL; do 
    # Loop through the different window sizes
    for d in 50 200 500; do 
        # Execute the R script with the specified parameters
        Rscript MAPS_create_control_file.R \
        -e /path/to/EEMS/${prefix}_eemsFiles/ \
        -i /path/to/MAPS/${prefix}_MAPSFiles \
        -p $prefix \
        -o TRUE \
        -n -1 \
        -d ${d} \
        -g 1000 \
        -l ${lwB} \
        -u ${upB} \
        -m 2000000 \
        -b 1000000 \
        -t 1000
    done 
done
