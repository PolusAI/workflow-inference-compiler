#!/bin/bash -e
cd ../..

if [ -d "image-workflows" ]
then
    echo "Directory image-workflows already exists."
    echo "Please install image-workflows manually."
    exit 1
else
    git clone https://github.com/PolusAI/image-workflows.git
fi
