#!/bin/bash
#echo "ROOT"
cat $1 | sed 's/\.$//g'
#echo "ENDROOT"
#echo "TORSDOF 0"