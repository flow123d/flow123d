#!/bin/bash
# compile_and_report.sh
# Usage: ./compile_and_report.sh [class_name] [source_file] [output_binary] [output_report] [compiler] [flags] [libs]

# Arguments
REPORT_TYPE=$1
CLASS_NAME=$2
OUTPUT_FILE=$3
shift 3
COMPILER_COMMAND=( "$@" )

LOG_TMP=$(mktemp)
SCRIPT_DIR="$(dirname "$(realpath "${BASH_SOURCE[0]}")")"

echo "Compiling ${CLASS_NAME} with type ${REPORT_TYPE} to generate report at ${OUTPUT_FILE}"
echo "Command: ${COMPILER_COMMAND}"

# Compile and capture output
"${COMPILER_COMMAND[@]}" > "$LOG_TMP" 2>&1
STATUS=$?


# Generate report
echo "Generating compilation report..."
python3 "$SCRIPT_DIR/runtime_and_compilation_reporter.py" \
    --status ${STATUS} \
    --output ${OUTPUT_FILE} \
    --log "$(cat ${LOG_TMP})" \
    --class-name ${CLASS_NAME} \
    --report-type ${REPORT_TYPE}

echo "Compilation status: ${STATUS}"
if [ ${STATUS} -eq 0 ]; then
    echo "Compilation of ${CLASS_NAME} successful. Report generated at ${OUTPUT_FILE}"
else
    echo "Compilation of ${CLASS_NAME} failed. Report generated at ${OUTPUT_FILE}"
fi

exit ${STATUS}