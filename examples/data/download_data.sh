#!/bin/bash

wget --no-clobber https://scs.illinois.edu/system/files/inline-files/NCIDiversity.rar
wget --no-clobber https://scs.illinois.edu/system/files/inline-files/NCINaturalChallenge.rar
wget --no-clobber https://scs.illinois.edu/system/files/inline-files/NCIOpen.rar

# bsdtar is installed with libarchive
bsdtar -xf NCIDiversity.rar
bsdtar -xf NCINaturalChallenge.rar
bsdtar -xf NCIOpen.rar
