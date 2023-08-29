#!/bin/bash -e
cd ../..

if [ -d "timeseriesplots" ]
then
    echo "Directory timeseriesplots already exists."
    echo "Please install timeseriesplots manually."
    exit 1
else
    git clone https://github.com/jfennick/timeseriesplots.git
    cd timeseriesplots
    cd install
    ./install_system_deps.sh
    cd ..
    pip install -e ".[all]"
fi
