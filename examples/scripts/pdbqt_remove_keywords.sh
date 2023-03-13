#!/bin/bash -e
grep -v TITLE "$1" | grep -v CRYST | grep -v REMARK | grep -v TER | grep -v MODEL | grep -v ENDMDL