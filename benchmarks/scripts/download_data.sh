#!/usr/bin/env bash
#
# This script downloads the required input data (bigWig tracks and BED peak
# sets) from the public S3 bucket where they are hosted.
#
# You must have the AWS CLI installed and configured for this to work.
#

set -euo pipefail

# Public S3 bucket URL where the data is hosted
BUCKET_URL="s3://iceqream"

# The directory where the data will be downloaded.
# Corresponds to the DATA_PREFIX variable in `commands.sh`.
DEST_DIR="benchmarks/data"

echo "––– Dataset Downloader –––"
echo "Downloading files from $BUCKET_URL to $DEST_DIR/"

# Create the destination directory if it doesn't exist.
mkdir -p "$DEST_DIR"

# Download specific files and directories needed for the benchmark
echo "Downloading Mouse Gastrulation Data..."
aws s3 sync "$BUCKET_URL/bw-cell-type-no-tss/" "$DEST_DIR/bw-cell-type-no-tss/"
aws s3 cp "$BUCKET_URL/gastru_peaks_no_tss.bed" "$DEST_DIR/gastru_peaks_no_tss.bed"
aws s3 sync "$BUCKET_URL/bw/" "$DEST_DIR/bw/"
aws s3 cp "$BUCKET_URL/gastru_peaks.bed" "$DEST_DIR/gastru_peaks.bed"

echo "Downloading Human Hematopoiesis Data..."
aws s3 sync "$BUCKET_URL/bw-bone-marrow/" "$DEST_DIR/bw-bone-marrow/"
aws s3 cp "$BUCKET_URL/bm_peaks.bed" "$DEST_DIR/bm_peaks.bed"

echo ""
echo "✅ All required benchmark data files have been downloaded."

 