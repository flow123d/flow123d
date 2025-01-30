#!/bin/bash

# Příkaz ke kompilaci jako jeden argument
COMPILE_CMD="$1"
COMPILE_OUTPUT="$2"      # Cesta k výstupnímu souboru kompilace
XML_REPORT="$3"          # Cesta k XML reportu


COMPILE_LOG="${COMPILE_OUTPUT}.log"

{
    eval "$COMPILE_CMD"
} &> "$COMPILE_LOG"

STATUS=$?

printf "Compile status: %d\n" "$STATUS"
python3 "$(dirname "$0")/generate_compile_xml.py" "$STATUS" "$COMPILE_LOG" "$XML_REPORT"

exit $STATUS
