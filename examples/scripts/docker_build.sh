#!/bin/bash -e
# docker does not follow symlinks, so copy
# See https://github.com/moby/moby/issues/1676
mkdir build/
cp "$(basename "$1")" "$(basename "$3")" build/
cd build/
# NOTE: basename will strip any username/ prefix, but this is intentional.
sudo docker build -f "$(basename "$1")" -t "$(basename "$2")" .

# NOTE: The nested double quotes above are apparently correct.
# See https://unix.stackexchange.com/questions/100945/word-splitting-when-parameter-is-used-within-command-substitution