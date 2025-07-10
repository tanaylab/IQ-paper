#!/usr/bin/env bash
#
# This script downloads the hg38 and mm10 reference genomes and chromosome
# size files from UCSC. It will place them in the current directory.
#
# The script checks if the files already exist before downloading,
# so it's safe to run multiple times.

set -euo pipefail

echo "––– Genome Downloader –––"
echo "This script will download required genome files into the current directory."
echo ""

# Human (hg38)
if [ ! -f "hg38.fa" ]; then
    echo "Downloading hg38 FASTA (hg38.fa)..."
    wget -q --show-progress -O hg38.fa.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
    gunzip hg38.fa.gz
else
    echo "hg38.fa already exists, skipping."
fi

if [ ! -f "hg38.chrom.sizes" ]; then
    echo "Downloading hg38 chrom.sizes..."
    wget -q --show-progress -O hg38.chrom.sizes "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"
else
    echo "hg38.chrom.sizes already exists, skipping."
fi

echo ""

# Mouse (mm10)
if [ ! -f "mm10.fa" ]; then
    echo "Downloading mm10 FASTA (mm10.fa)..."
    wget -q --show-progress -O mm10.fa.gz http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
    gunzip mm10.fa.gz
else
    echo "mm10.fa already exists, skipping."
fi

if [ ! -f "mm10.chrom.sizes" ]; then
    echo "Downloading mm10 chrom.sizes..."
    wget -q --show-progress -O mm10.chrom.sizes "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes"
else
    echo "mm10.chrom.sizes already exists, skipping."
fi

echo ""
echo "✅ All genome files are present."
 