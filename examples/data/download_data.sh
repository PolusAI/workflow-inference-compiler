#!/bin/bash

# NCI databases
wget --no-clobber https://scs.illinois.edu/system/files/inline-files/NCIDiversity.rar
wget --no-clobber https://scs.illinois.edu/system/files/inline-files/NCINaturalChallenge.rar
wget --no-clobber https://scs.illinois.edu/system/files/inline-files/NCIOpen.rar

# bsdtar is installed with libarchive
bsdtar -xf NCIDiversity.rar
bsdtar -xf NCINaturalChallenge.rar
bsdtar -xf NCIOpen.rar

# SMACC database
wget --no-clobber https://smacc.mml.unc.edu/ncats_phenotypic_curated.xlsx
wget --no-clobber https://smacc.mml.unc.edu/ncats_target_based_curated.xlsx
