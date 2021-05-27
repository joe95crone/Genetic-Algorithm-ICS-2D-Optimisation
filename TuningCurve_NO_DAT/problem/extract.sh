#!/bin/sh

rms=`grep Bandwidth ColFlux_Bandwidth_Calc.txt |awk '{print $3}'`
fl=`grep Collimated ColFlux_Bandwidth_Calc.txt |awk '{print $3}'`
echo $rms -$fl
