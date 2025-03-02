#!/bin/bash
# Usage: ./compile_and_report.sh <xml_output> <compile command and his arguments>

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <xml_output> <compilation_command>"
    exit 1
fi

echo "Xml output: $1"
echo "Compilation command: $2"

XML_OUTPUT="$1"
COMMANDS="$@"

LOG_TMP=$(mktemp)

COMMANDS > "$LOG_TMP" 2>&1
STATUS=$?

python3 compile_reporter.py --status "$STATUS" --output "$XML_OUTPUT" --log "$(cat "$LOG_TMP")"

rm "$LOG_TMP"
exit "$STATUS"
