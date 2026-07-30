#!/bin/bash

# Check if a version number is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: bash $0 <version_number>"
    exit 1
fi

# Version number from the first argument
VERSION=$1

# Directory for the release
RELEASE_DIR="../releases"
VERSION_DIR="$RELEASE_DIR/polegon-$VERSION-alpha-linux-x86_64"

# Create version directory
mkdir -p $VERSION_DIR

# Compile the program with optimizations and debugging information
g++ -std=c++17 -O3 -g -static -fopenmp *.cpp -o $VERSION_DIR/polegon

# Compile the debug version of the program
g++ -std=c++17 -g -static -fopenmp *.cpp -o $VERSION_DIR/polegon_debug

# Copy additional files
cp $VERSION_DIR/polegon polegon
cp $VERSION_DIR/polegon_debug polegon_debug
cp polegon_master $VERSION_DIR/polegon_master
cp ../LICENSE $VERSION_DIR/LICENSE

# Change directory to releases
cd $RELEASE_DIR

# Create a tarball with the version number in the name
tar -cvzf "polegon-$VERSION-alpha-linux-x86_64.tar.gz" "polegon-$VERSION-alpha-linux-x86_64" 
rm -rf "polegon-$VERSION-alpha-linux-x86_64" 
