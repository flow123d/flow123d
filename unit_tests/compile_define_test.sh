#!/bin/bash

echo "DEBUG: Script started"
echo "Compile command: $1"
echo "Compile output: $2"
echo "XML report: $3"

COMPILE_CMD="$1"         # Commands to compile
COMPILE_OUTPUT="$2"      # Path to output file
XML_REPORT="$3"          # Path to XML report

COMPILE_CMD=$(echo "$COMPILE_CMD" | sed 's/;/ /g')

{
    eval "$COMPILE_CMD"
} > "$COMPILE_OUTPUT" 2>&1

STATUS=$?
printf "Compile status: %d\n" "$STATUS"

BINARY_PATH=$(echo "$COMPILE_CMD" | awk '{for (i=1; i<=NF; i++) if ($i == "-o") print $(i+1)}')

if [ -f "$BINARY_PATH" ]; then
    echo "Binary exists: $BINARY_PATH"
    ls -l "$BINARY_PATH"
else
    echo "Binary not found: $BINARY_PATH"
    STATUS=127
fi

echo "DEBUG: Starting python script"

python3 "$(dirname "$0")/generate_compile_xml.py" "$STATUS" "$COMPILE_OUTPUT" "$XML_REPORT"

echo "DEBUG: Finished python script"

exit $STATUS


