#!/bin/bash
# Usage: ./compile_and_report.sh <xml_report_file> <compilation command and its arguments>

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <xml_report_file> <compilation_command>"
    exit 1
fi

XML_REPORT="$1"
shift
COMMANDS="$@"
LOG_TMP=$(mktemp)

echo "XML report will be generated at: ${XML_REPORT}"
echo "Compilation command: $@"

eval "$COMMANDS" > "$LOG_TMP" 2>&1
STATUS=$?

python3 "$(dirname "$0")/generate_jtest_xml.py" --status "$STATUS" --output "$XML_REPORT" --log "$(cat "$LOG_TMP")"

rm "$LOG_TMP"
exit "$STATUS"
