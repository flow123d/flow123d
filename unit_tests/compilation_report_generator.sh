#!/bin/bash

# Script to generate compilation report XML
# Usage: compile_report_generator.sh [class_name] [log_file] [output_xml]

CLASS_NAME=$1
LOG_FILE=$2
OUTPUT_XML=$3

# Check if the binary exists - if it does, build was successful
TEST_BINARY="${CMAKE_CURRENT_BINARY_DIR}/${CLASS_NAME}_test_bin"
if [ -f "$TEST_BINARY" ] || [ -f "${TEST_BINARY}.exe" ]; then
    BUILD_STATUS=0
else
    BUILD_STATUS=1
fi

# Get the log content
LOG_CONTENT=$(cat "$LOG_FILE")

# Call the Python script to generate the report
python3 compilation_reporter.py --status $BUILD_STATUS --output "$OUTPUT_XML" --log "$LOG_CONTENT" --class "$CLASS_NAME"

# Print status message
if [ $BUILD_STATUS -eq 0 ]; then
    echo "Compilation of $CLASS_NAME successful. Report generated at $OUTPUT_XML"
else
    echo "Compilation of $CLASS_NAME failed. Report generated at $OUTPUT_XML"
fi

exit $BUILD_STATUS