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
    ./conda_devtools.sh
    pip install -e ".[all]"
fi

cd ../workflow-inference-compiler/install/
