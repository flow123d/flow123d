#!/bin/bash

echo "DEBUG: Script started"
echo "Compile command: $1"
echo "Compile output: $2"
echo "XML report: $3"

COMPILE_CMD="$1"         # Commands to compile
COMPILE_OUTPUT="$2"      # Path to output file
XML_REPORT="$3"          # Path to XML report

mkdir -p "$(dirname "$XML_REPORT")"

{
    eval "$COMPILE_CMD"
} > "$COMPILE_OUTPUT" 2>&1

STATUS=$?
printf "Compile status: %d\n" "$STATUS"

# If compilation fails, we generate an XML report
if [ $STATUS -ne 0 ]; then
    python3 "$(dirname "$0")/generate_compile_xml.py" "$STATUS" "$COMPILE_OUTPUT" "$XML_REPORT"
    exit $STATUS
fi

# If compilation is successful, but the binary is missing, we also report it
if [ ! -f "$(dirname "$COMPILE_OUTPUT")/$(basename "$COMPILE_OUTPUT" _compile_output.txt)" ]; then
    echo "Error: Test binary missing!" >> "$COMPILE_OUTPUT"
    python3 "$(dirname "$0")/generate_compile_xml.py" 127 "$COMPILE_OUTPUT" "$XML_REPORT"
    exit 127
fi

exit $STATUS

