#!/bin/bash
# compile_and_report.sh
# Usage: ./compile_and_report.sh <class_name> <output_file> <compiler_command>

CLASS_NAME=$1
OUTPUT_FILE=$2
shift 2
COMPILER_COMMAND=("$@")

LOG_TMP=$(mktemp)
SCRIPT_DIR="$(dirname "$(realpath "${BASH_SOURCE[0]}")")"

echo "Compiling ${CLASS_NAME} to run tests and generate report at ${OUTPUT_FILE}"
echo "Command: ${COMPILER_COMMAND[@]}"

"${COMPILER_COMMAND[@]}" > "$LOG_TMP" 2>&1
STATUS=$?

echo "Compilation status: ${STATUS}"
if [ ${STATUS} -eq 0 ]; then
    echo "Compilation of ${CLASS_NAME} successful."
else
    echo "Compilation of ${CLASS_NAME} failed. Report generated at ${OUTPUT_FILE}"
    echo "Generating compilation report..."
    python3 "$SCRIPT_DIR/runtime_and_compilation_reporter.py" \
        --output ${OUTPUT_FILE} \
        --log "$(cat ${LOG_TMP})" \
        --class-name ${CLASS_NAME}
fi

exit ${STATUS}