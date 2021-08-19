#!/bin/sh

pf=`grep Penalty ColFlux_penfunc_Calc.txt |awk '{print $3}'`
fl=`grep Collimated ColFlux_penfunc_Calc.txt |awk '{print $3}'`
echo $pf -$fl
