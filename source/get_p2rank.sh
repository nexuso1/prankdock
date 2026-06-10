#!/bin/bash

# Script to download and extract p2rank

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Parent directory (the repo root)
REPO_DIR="$(dirname "$SCRIPT_DIR")"
# Target directory for p2rank
P2RANK_DIR="$REPO_DIR/p2rank"

# Create the p2rank directory if it doesn't exist
mkdir -p "$P2RANK_DIR"

# Download the p2rank archive
echo "Downloading p2rank 2.5.1..."
wget -O /tmp/p2rank_2.5.1.tar.gz https://github.com/rdk/p2rank/releases/download/2.5.1/p2rank_2.5.1.tar.gz

# Check if download was successful
if [ $? -ne 0 ]; then
    echo "Error: Failed to download p2rank"
    exit 1
fi

# Extract the archive
echo "Extracting p2rank..."
tar -xzf /tmp/p2rank_2.5.1.tar.gz -C "$P2RANK_DIR" --strip-components=1

# Check if extraction was successful
if [ $? -ne 0 ]; then
    echo "Error: Failed to extract p2rank"
    exit 1
fi

# Clean up the downloaded archive
rm /tmp/p2rank_2.5.1.tar.gz

echo "p2rank has been successfully downloaded and extracted to $P2RANK_DIR"
