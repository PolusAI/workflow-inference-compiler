#!/bin/bash
echo "ROOT"
cat $1 | grep -v ROOT | grep -v BRANCH | grep -v TORSDOF
echo "ENDROOT"
echo "TORSDOF 0"