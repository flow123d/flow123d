# Configuration for CI server
# debug build

set(FLOW_BUILD_TYPE debug)

set(CMAKE_VERBOSE_MAKEFILE on)

set(PETSC_INSTALL_CONFIG "bddcml")
set(BDDCML_ROOT "_INSTALL_")
set(USE_PYTHON "yes")

set(PLATFORM_NAME "linux_x86_64")


# --------------------------------------
set(PYTHON_COPY_ROOT /usr/lib64/python2.7)
set(PYTHON_COPY_PACKAGES
    bsddb compiler ctypes curses distutils email encodings hotshot idlelib 
    importlib json lib-dynload lib-tk lib2to3 logging multiprocessing 
    plat-linux2 pydoc_data sqlite3 test unittest wsgiref xml)
set(PYTHON_SYSPATH
    plat-linux2 lib-tk lib-old lib-dynload dist-packages site-packages)
set(PYTHON_COPY_DIST_PACKAGES)