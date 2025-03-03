#!/bin/bash
# compile_and_report.sh
# Usage: ./compile_and_report.sh [class_name] [source_file] [output_binary] [output_report] [compiler] [flags] [libs]

# Arguments
CLASS_NAME=$1
SOURCE_FILE=$2
OUTPUT_BINARY=$3
OUTPUT_REPORT=$4
COMPILER=$5
FLAGS=$6
LIBS=$7

# Define log files
LOG_OUTPUT="${OUTPUT_BINARY}_compile_output.log"
LOG_ERRORS="${OUTPUT_BINARY}_compile_errors.log"
LOG_COMBINED="${OUTPUT_BINARY}_combined.log"

echo "Compiling ${CLASS_NAME} test..."
echo "Command: ${COMPILER} ${FLAGS} -o ${OUTPUT_BINARY} ${SOURCE_FILE} ${LIBS}"

# Compile and capture output
${COMPILER} ${FLAGS} -o ${OUTPUT_BINARY} ${SOURCE_FILE} ${LIBS} > ${LOG_OUTPUT} 2> ${LOG_ERRORS}
COMPILE_STATUS=$?

# Combine logs
cat ${LOG_OUTPUT} ${LOG_ERRORS} > ${LOG_COMBINED}

# Generate report
echo "Generating compilation report..."
python3 "$(pwd)/compilation_reporter.py" \
    --status ${COMPILE_STATUS} \
    --output ${OUTPUT_REPORT} \
    --log "$(cat ${LOG_COMBINED})"

echo "Compilation status: ${COMPILE_STATUS}"
if [ ${COMPILE_STATUS} -eq 0 ]; then
    echo "Compilation of ${CLASS_NAME} successful. Report generated at ${OUTPUT_REPORT}"
else
    echo "Compilation of ${CLASS_NAME} failed. Report generated at ${OUTPUT_REPORT}"
fi

exit ${COMPILE_STATUS}