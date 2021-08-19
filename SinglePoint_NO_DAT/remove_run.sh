#!/bin/bash

# remove all extraneous files for new GA run
rm *png # removes GA plots
(cd problem; rm -r run*) # removes GA problem directories
(cd run; rm *out *txt *log)


