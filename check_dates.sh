#!/bin/bash
WORK_DIR="/home/jb/workspace/flow123d"
DATE_FILE="${WORK_DIR}/check_date"
if [ -f "${DATE_FILE}" ]
then
	find "${WORK_DIR}" -newer "${DATE_FILE}" -type f -print
fi
touch "${DATE_FILE}"

