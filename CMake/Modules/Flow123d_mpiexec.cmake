####
# this CMake script is called at build time to make calling script to given mpiexec binary
# assume variable PETSC_MPIEXEC, 


if(EXISTS ${PETSC_MPIEXEC})
  file(REMOVE ${CMAKE_BINARY_DIR}/bin/mpiexec)
  file(WRITE ${CMAKE_BINARY_DIR}/bin/mpiexec "#!/bin/bash\n")
  file(APPEND ${CMAKE_BINARY_DIR}/bin/mpiexec "\"${PETSC_MPIEXEC}\" \"\$@\"")
else()
  if(COMMAND mpiexec)
    message(WARNING "Missing mpiexec in PETSc instalation. Using system wide mpiexec.")
    file(REMOVE ${CMAKE_BINARY_DIR}/bin/mpiexec)
    file(WRITE ${CMAKE_BINARY_DIR}/bin/mpiexec "#!/bin/bash\n")
    file(APPEND ${CMAKE_BINARY_DIR}/bin/mpiexec "mpiexec \"\$@\"")
  else()
    message(WARNING "Missing any mpiexec.")
  endif()
endif()

