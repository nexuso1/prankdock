#!/bin/bash

# Script to download and extract p2rank

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Parent directory (the repo root)
REPO_DIR="$(dirname "$SCRIPT_DIR")"
# Target directory for p2rank
P2RANK_DIR="$REPO_DIR/p2rank"
CAVER_DIR="$REPO_DIR/caver_3.0"
CAVERDOCK_PATH="$REPO_DIR/caverdock-1.2.sif"

mkdir -p "$REPO_DIR"/tmp

# Create the p2rank directory if it doesn't exist
mkdir -p "$P2RANK_DIR"

# Download the p2rank archive
echo "Downloading p2rank 2.5.1..."
wget -O -q --show_progress "$REPO_DIR"/tmp/p2rank_2.5.1.tar.gz https://github.com/rdk/p2rank/releases/download/2.5.1/p2rank_2.5.1.tar.gz

# Check if download was successful
if [ $? -ne 0 ]; then
    echo "Error: Failed to download p2rank"
    exit 1
fi

# Extract the archive
echo "Extracting p2rank..."
tar -xzf "$REPO_DIR"/tmp/p2rank_2.5.1.tar.gz -C "$P2RANK_DIR" --strip-components=1

# Check if extraction was successful
if [ $? -ne 0 ]; then
    echo "Error: Failed to extract p2rank"
    exit 1
fi

# Clean up the downloaded archive
rm /tmp/p2rank_2.5.1.tar.gz

echo "p2rank has been successfully downloaded and extracted to $P2RANK_DIR"

# Download CaverDock
echo "Downloading CaverDock ..."
wget -O -q --show-progress "$CAVERDOCK_PATH" https://loschmidt.chemi.muni.cz/static/releases/caverdock/1.2/caverdock-1.2.sif
if [ $? -ne 0 ]; then
    echo "Error: Failed to download CaverDock"
    exit 1
fi
echo "Done. Saved to $CAVERDOCK_PATH"

# Download CAVER
echo "Downloading CAVER..."
wget -O -q --show-progress "$REPO_DIR"/tmp/caver_3.0.zip https://www.caver.cz/fil/download/caver30/301/caver_3.0.zip
echo "Extracting CAVER..."
mkdir -p "$CAVER_DIR"
unzip "$REPO_DIR"/tmp/caver_3.0.zip -d "$CAVER_DIR"

if [ $? -ne 0 ]; then
    echo "Error: Failed to extract CAVER"
    exit 1
fi

rm "$REPO_DIR"/tmp/caver_3.0.zip

echo "Done. Extracted to $CAVER_DIR"

conda create --file environment.yaml