INI_FILES="./flow.ini ./flow_vtk.ini"
NPROC="1"
FLOW_PARAMS="-i input -o output -ksp_atol 1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor"

test:
	../../bin/tests/run_test.sh ${INI_FILES} ${NPROC} ${FLOW_PARAMS}

clean:
	rm -f profiler_*
	rm -f output/*