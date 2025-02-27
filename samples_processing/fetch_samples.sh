#!/bin/bash

# Define variables
DATASET_ID="</projects/proj-6030_ntrk2isoforms-1128.4.1092/samples>"  # Mediaflux path
SAMPLES=("10X04_4.bam.1" "10X52_2.bam.1" "10X26_4.bam.1" "10X04_2.bam.1" "10X26_3.bam.1")
OUTPUT_DIR="samples"

# Create output directory
mkdir -p $OUTPUT_DIR

# Download each sample
for SAMPLE in "${SAMPLES[@]}"; do
  mflux dataset fetch $DATASET_ID/$SAMPLE --output $OUTPUT_DIR/$SAMPLE
done
