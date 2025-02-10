#!/bin/bash

# Příkaz ke kompilaci jako jeden argument
COMPILE_CMD="$1"
COMPILE_OUTPUT="$2"      # Cesta k výstupnímu souboru kompilace
XML_REPORT="$3"          # Cesta k XML reportu


{
    eval "$COMPILE_CMD"
} &> "$COMPILE_OUTPUT"

STATUS=$?

printf "Compile status: %d\n" "$STATUS"
python3 "$(dirname "$0")/generate_compile_xml.py" "$STATUS" "$COMPILE_OUTPUT" "$XML_REPORT"

exit $STATUS
