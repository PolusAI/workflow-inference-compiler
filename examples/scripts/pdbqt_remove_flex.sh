#!/bin/bash -e
echo "ROOT"
grep -v ROOT "$1" | grep -v BRANCH | grep -v TORSDOF
echo "ENDROOT"
echo "TORSDOF 0"