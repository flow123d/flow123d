#!/bin/bash

set -x
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
FLOW_ROOT_DIR="${SCRIPT_DIR}/.."

# Script to run a benchmark unit test with defined number of repetitions

function print_usage() {
cat << EOF
  
Usage:  
  
    run_benchmark.sh [-np=n_proc] [-nr=n_runs] [-t=timeout] [--help] (executable|test_path)

    
n_runs could be provided by the system variable BENCHMARK_N_RUNS, comand line argument takes the priority

Test could be specified either as <executable>, i.e. path to the test binary, or as a <test_path> in form:
"<unit_test_dir>/<test_target_base name>". For example "coupling/dg_asm".
    
EOF

}

function arg_assignment_split() {
    RESULT_ARG=${1%=*}
    if [ "$RESULT_ARG" == "$1" ]
    then  
        RESULT_VALUE=
    else
        RESULT_VALUE="${1#*=}"
    fi
}

function parse_arguments() {
    N_RUNS=${BENCHMARK_N_RUNS:-1}
    N_PROC=1
    TIMEOUT=300      # default timeout 5 min
    
    while [ "${1#-}" != "$1" ]      # arg starts with '-'
    do
    arg_assignment_split "$1"       # produce $RESULT_ARG and $RESULT_VALUE
    arg="$RESULT_ARG"
    value="$RESULT_VALUE"
    shift
    case $arg in
        -np)
        N_PROC=${value}
        ;;  
        -nr)        
        if [ -z "${value}" ]; then
            echo "-nr missing value"
            exit 1
        fi
        N_RUNS=${value}
        ;;
        -t)
        TIMEOUT=${value}
        ;;
        -h|--help)
        print_usage
        exit 0
        ;;
        *)
        print_usage
        error "Invalid argument '$arg'"
        ;;
    esac
    done
    
    test_spec="$1"
    shift
}

parse_arguments "$@"

GIT_SHORT_HASH=`git rev-parse --short HEAD`

BUILD_DIR="${FLOW_ROOT_DIR}/build_tree"

if [ -x "${test_spec}" ]; then
    test_binary="${test_spec}"
    file_name="${test_binary##*/}"
    dir="${test_binary%/*}"
    TEST_TARGET_BASE=${file_name%_bench_bin}
    TEST_SUBDIR=${dir##*/}
else
    # test_spec is a relative test target path in format "<unit_test_dir>/<test_target_base name>"
    TEST_SUBDIR="${test_spec%/*}"
    TEST_TARGET_BASE="${test_spec#*/}"
fi

TEST_ABS_DIR="${BUILD_DIR}/unit_tests/${TEST_SUBDIR}"
TEST_BINARY=${TEST_ABS_DIR}/${TEST_TARGET_BASE}_bench_bin
if [ ! -x ${TEST_BINARY} ];then
    echo "test_spec=${test_spec} leads to unknown test_binary: ${TEST_BINARY}"
    exit 1
fi

cd ${TEST_ABS_DIR}
for (( run_id=1; run_id<=$N_RUNS; run_id++ ));do
    TEST_NAME=${TEST_TARGET_BASE}-${N_PROC}-bench
    ${SCRIPT_DIR}/time_limit.sh -t ${TIMEOUT} ${BUILD_DIR}/bin/mpiexec -np ${N_PROC} ${TEST_BINARY} --gtest_output="xml:${TEST_NAME}.xml"
    ls -l
    mv ${TEST_TARGET_BASE}_profiler.json ${TEST_NAME}_profiler_${run_id}.json
    # python3 ${FLOW_ROOT_DIR}/unit_tests/${TEST_SUBDIR}/${TEST_TARGET_BASE}_bench_postprocess.py ${TEST_NAME}_profiler.json ${GIT_SHORT_HASH} ${run_id}
done
