#!/bin/bash

gen=`grep GENERATION run/var_output.txt|tail -1`
best=`sort -u -k 2 run/var_output.txt|head -2|tail -1`

best2=`grep -B 100 "GENERATION $1" run/var_output.txt|sort -u -k 2| head -2 |tail -1`
echo "Best up to " $gen $best
echo "Best up to  ### GENERATION" $1 $best2
