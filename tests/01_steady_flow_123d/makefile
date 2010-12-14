INI_FILES="/cygdrive/c/cygwin_drive/new/tests/01_steady_flow_123d/flow_schur2.ini"
NPROC="2"
FLOW_PARAMS=" -ksp_atol 1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor"

test:
	@cd ..;\
	bash ./run_test.sh ${INI_FILES} ${NPROC}${FLOW_PARAMS}
