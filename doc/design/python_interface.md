# Notes to creating Python interface for Flow123d

## Research 
Resources for creating the Python interface to a library:
[Static binary wheel approach "YOLO"](https://thomastrapp.com/posts/building-a-pypi-package-for-a-modern-cpp-project//)

Linking to libpython is [not recomended](https://www.python.org/dev/peps/pep-0513/#libpythonx-y-so-1)
setting RTLD_GLOBAL flag for dlopen would be necessary.

pybind11 - propsed setup.py includes building of the extension *.c sources

## Possible schema

- do not try to make python PIP package for Flow123d
- provide some installation process
- use it to install Flow123d and python interface in distributed images

Install process:
- use pybind11 example setup.py with cmake support
- try to separate build of flow123d libs from build of the Python interface
- call the build of python interface from setup.py
