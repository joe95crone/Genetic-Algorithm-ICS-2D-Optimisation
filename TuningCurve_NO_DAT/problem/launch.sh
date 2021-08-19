#!/bin/bash

sed "s%<THETA>%$1%g;s%<BETAX>%$2%g;s%<BETAY>%$3%g;" ../template/indata.temp > indata.txt
nice python3 ../../ICS_FBW.py
