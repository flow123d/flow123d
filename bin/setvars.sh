#!/bin/bash
FILE=/opt/intel/oneapi/setvars.sh
if [ -f "$FILE" ]; then
    . $FILE
fi
