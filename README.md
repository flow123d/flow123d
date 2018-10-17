---
title: Overview
layout: page
---


# Flow123d
*Transport Processes in Fractured Media*

<!-- build status -->
[![Jenkins release](https://img.shields.io/jenkins/s/http/ciflow.nti.tul.cz:8080/Flow123d-ci2runner-release-multijob.svg?style=flat-square&label=jenkins%20release)](http://ciflow.nti.tul.cz:8080/view/multijob-list/job/Flow123d-ci2runner-release-multijob/)
[![Jenkins debug](https://img.shields.io/jenkins/s/http/ciflow.nti.tul.cz:8080/Flow123d-ci2runner-debug-multijob.svg?style=flat-square&label=jenkins%20debug)](http://ciflow.nti.tul.cz:8080/view/multijob-list/job/Flow123d-ci2runner-debug-multijob/)
[![Coveralls master](https://img.shields.io/coveralls/github/flow123d/flow123d.svg?style=flat-square&label=coverage%20master)](https://coveralls.io/github/flow123d/flow123d)

Flow123d is a simulator of underground water flow and transport processes in fractured
porous media. Novelty of this software is support of computations on complex
meshes consisting of simplicial elements of different dimensions. Therefore
we can combine continuum models and discrete fracture network models.
For more information see the project pages:
[flow123d.github.com](http://flow123d.github.com).


## Installation
This guide summarises how to install Flow123d. The easiest and the most
reliable way to get started is to use Docker container (available for Linux and Windows).

For installation in high performance computing clusters we recommend always building **from source**.

*Note: Windows platform support is limited and is therefore recommended to use
stable only packages for the installation.*


### From source
Flow123d requires several libraries and tools.
The following can be easily installed via packaging tools (such as `apt` or `yum`):
 - [CMake](https://cmake.org/), make, gcc, g++, python3, perl, gfortran, git and Boost.

Other libraries can be a bit trickier to install:
 - [Yamlcpp](https://github.com/jbeder/yaml-cpp),
   [Armadillo](http://arma.sourceforge.net/),
   [PETSc](https://www.mcs.anl.gov/petsc/),
   [BDDCML](http://users.math.cas.cz/sistek/software/bddcml.html) and
   [Blopex](https://bitbucket.org/joseroman/blopex).

To see what libraries are precisely required, head over to
[flow123d-docker-images](https://github.com/flow123d/flow123d-docker-images) repository
where you can learn more.

---

In order to install Flow123d from the source, first, clone the repository and enter the
main directory:
```sh
git clone https://github.com/flow123d/flow123d.git
cd flow123d/
```
Create a `config.cmake` configuration file (see [config.cmake.template](config.cmake.template) file).
For an inspiration look at various `config.cmake` files in the [/config](/config) folder.
To start a compilation, simply run:
```sh
make
```
To speed up the entire compilation process, add flag `-j <cores>`, where `<cores>` is
a number of cores, you want to use. For example if your PC has 16 cores, run the following:
```sh
make -j 16
```

### Containers/Docker
#### Linux

To get started, install [Docker](https://www.docker.com/) and clone the
repository, then enter the main directory and copy the docker configuration file:
```sh
git clone https://github.com/flow123d/flow123d.git
cd flow123d/
cp /config/config-jenkins-docker-debug.cmake config.cmake
```
*Note: to use release version of the libraries, substitute `debug` for the `release`*

After the configuration has been prepared, enter the docker image via utility wrapper
script `fterm`:

> An `fterm` script is a wrapper for the docker. It is recommended to
> use this tool instead of raw docker calls as it automatically mounts necessary paths
> and takes care of various other necessities such as perserving file permission.

```sh
bin/fterm
```

*Note: by default, `fterm` is using a docker image which contains a `debug` version of libraries.
They can be useful to a developer but for production use you may want to switch to a `release` version.
To do that, first copy a `config-jenkins-docker-release.cmake` file instead of a `debug` one.
Next, run `fterm` with a flag `rel` (as for a `release`):*

```sh
bin/fterm rel
```

This process may take **several minutes** for the first time
(it greatly depends on your internet connection). It will download a debug (or a release)
docker image from the [docker hub](https://hub.docker.com/u/flow123d/) and personalize
the image for the current user (This will preserve files permissions and won't cause
any troubles with root access combined with MPI). To see more of the `fterm`'s capabilities
execute `bin/fterm --help`.

After the command, you should enter the docker container and see something like this:
```sh
 ___ _            _ ___ ____    _  
| __| |_____ __ _/ |_  )__ / __| |
| _|| / _ \ V  V / |/ / |_ \/ _  |
|_| |_\___/\_/\_/|_/___|___/\__,_|
                         dbg-3.0.0
jan-hybs@dbg-3.0.0 ~/projects/flow123d/flow123d
```

You should be now all set up for the compilation, simply run `make` while inside
the docker container and wait for the compilation to finish. This may
take up to **15 minutes**, depending on your computer's performance. You can
speed up the compilation by specifying the `-j` flag
(See section [From source](#from-source))


#### Windows
To use Flow123d under Windows, you have to first install [Docker Toolbox](https://www.docker.com/products/docker-toolbox).

*Note: If Docker Toolbox is not installed it can be automatically installed by the
Windows Flow123d installator, but the installator may not be
packing the latest version fo the Docker Toolbox.*

The installator uses `powershell` to start installtion process, make sure the binary is in the system path.
Download the installator from the official [website](http://flow123d.github.io/) and unpack it. Double click
`install.bat` and follow the instructions on the screen. *For the first time, you will be prompted
to allow script to execute.*


## Getting started
Please refer to a **User Guide and Input Reference manual** available
at our [official website](http://flow123d.github.io/) where there is a entire section dedicated
to this topic. You can find step-by-step tutorial explaining geometries, `yaml` input files
and more. Below you can see a result from the tutorial problem.
![](/doc/graphics/figure.png)



## Developers
### Troubleshooting
When problem occures during the compilation process it may be due to a leftover files in a build folder.
Cleaning this directory can solve this issue. You can either remove `build-<branch>` folder
(the folder is located one level above repository root) via
`make clean-all`, which removes build folders and also remove any symlinks.

To clean all the build fodlers manually run `rm -rf ../build-*` while in a repository root.
**Running `rm -rf` can quite easily cause a lot of damage, double check that you're
in a correct folder.**



### Building the reference manual

The reference manual can be built by while in docker container
```sh
make ref-doc
```
To copy out reference manual from docker use command
[`docker cp`](https://docs.docker.com/engine/reference/commandline/cp/).
