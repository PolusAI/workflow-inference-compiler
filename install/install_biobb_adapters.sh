#!/bin/bash -e
cd ../..

if [ -d "biobb_adapters" ]
then
    echo "Directory biobb_adapters already exists."
    echo "Please install biobb_adapters manually."
    exit 1
else
    # NOTE: Do not use https://github.com/bioexcel/biobb_adapters
    # Our fork has a few changes from the upstream repo
    git clone https://github.com/sameeul/biobb_adapters.git
fi
