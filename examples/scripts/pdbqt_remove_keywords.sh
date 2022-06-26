#!/bin/bash
cat $1 | grep -v TITLE | grep -v CRYST | grep -v REMARK | grep -v TER | grep -v MODEL | grep -v ENDMDL