#!/bin/bash
# Simple script which will use correct python version and execute runtest.py

# populated by configure_file call
/usr/bin/python2.7 /home/jb/workspace/flow123d/bin/python/runtest.py $@
