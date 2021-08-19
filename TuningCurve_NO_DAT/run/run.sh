#!/bin/bash

nice ../var/var var.conf ../pisa/ 0.2 < /dev/null > var.out 2>&1 &
sleep 5
nice ../spea2/spea2 spea2.conf ../pisa/ 0.2 < /dev/null > spea2.out 2>&1 &
