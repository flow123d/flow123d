INI_FILES="flow_schur2.ini"
NPROC="1 2"
TEST_DIR="01_steady_flow_123d"
FLOW_PARAMS=" -ksp_atol 1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor"

test:
	@cd ..;\
	bash ./run_test.sh ${INI_FILES} ${NPROC} ${TEST_DIR} ${FLOW_PARAMS}
