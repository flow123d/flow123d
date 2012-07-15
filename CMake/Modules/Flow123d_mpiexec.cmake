####
# this CMake script is called at build time to make calling script to given mpiexec binary
# assume variable PETSC_MPIEXEC, 


if(EXISTS ${PETSC_MPIEXEC})
  set(MPIEXEC_PATH ${PETSC_MPIEXEC})
else()
  if(COMMAND mpiexec)
    message(WARNING "Missing mpiexec in PETSc instalation. Using system wide mpiexec.")
    set(MPIEXEC_PATH mpiexec)
  else()
    message(WARNING "Missing any mpiexec.")
  endif()
endif()

configure_file(${CMAKE_SOURCE_DIR}/CMake/mpiexec_link_template ${CMAKE_BINARY_DIR}/mpiexec)

