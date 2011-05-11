INI_FILES="./flow.ini ./flow_vtk.ini"
#INI_FILES="./flow.ini"
NPROC="1 2 3"
FLOW_PARAMS=" -ksp_atol 1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor -i input -o output"

test:
	bash ../run_test.sh ${INI_FILES} ${NPROC} ${FLOW_PARAMS}
