#!/bin/bash

echo "DEBUG: Script started"
echo "Compile command: $1"
echo "Compile output: $2"
echo "XML report: $3"

COMPILE_CMD="$1"         # Commands to compile
COMPILE_OUTPUT="$2"      # Path to output file
XML_REPORT="$3"          # Path to XML report

COMPILE_CMD=$(echo "$COMPILE_CMD" | sed 's/;/ /g')
mkdir -p "$(dirname "$XML_REPORT")"

{
    eval "$COMPILE_CMD"
} 2>&1 | tee "$COMPILE_OUTPUT"

STATUS=$?

printf "Compile status: %d\n" "$STATUS"
cat "$COMPILE_OUTPUT"
python3 "$(dirname "$0")/generate_compile_xml.py" "$STATUS" "$COMPILE_OUTPUT" "$XML_REPORT"

exit $STATUS


STATUS=$?
