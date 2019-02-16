# Flow123d &middot; [![Jenkins release](https://img.shields.io/jenkins/s/http/ciflow.nti.tul.cz:8080/Flow123d-ci2runner-release-multijob.svg?style=flat-square&label=release)](http://ciflow.nti.tul.cz:8080/view/multijob-list/job/Flow123d-ci2runner-release-multijob/) [![Jenkins debug](https://img.shields.io/jenkins/s/http/ciflow.nti.tul.cz:8080/Flow123d-ci2runner-debug-multijob.svg?style=flat-square&label=debug)](http://ciflow.nti.tul.cz:8080/view/multijob-list/job/Flow123d-ci2runner-debug-multijob/) [![Coveralls master](https://img.shields.io/coveralls/github/flow123d/flow123d.svg?style=flat-square&label=coverage)](https://coveralls.io/github/flow123d/flow123d) [![Docker hub](https://img.shields.io/badge/docker-hub-blue.svg?colorA=2271b8&colorB=dc750d&logo=docker&style=flat-square&logoColor=FFF)](https://hub.docker.com/u/flow123d/) [![CI-HPC](https://img.shields.io/badge/ci--hpc-performace-green.svg?style=flat-square)](http://hybs.nti.tul.cz/ci-hpc/)

*Transport Processes in Fractured Media*


Flow123d is a simulator of underground water flow and transport processes in fractured
porous media. Novelty of this software is support of computations on complex
meshes consisting of simplicial elements of different dimensions. Therefore
we can combine continuum models and discrete fracture network models.
For more information see the project pages:
[flow123d.github.com](http://flow123d.github.com).


## Installation
This guide summarises how to install Flow123d. They are two ways to use Flow123d:

  1) [**Docker containers**](#docker-containers):
   - *This is the easiest and the most reliable way to start using Flow123d.
   It is available for Linux and Windows users.*  
   *Windows platform support is limited and is therefore recommended to use
   stable only packages for the installation.*

  2) [**From source**](#from-source):
   - *This option is only available to Linux users.
   It requires a moderate to high understanding of a Linux operating system.*


**For installation in high performance computing clusters we recommend always building from source**.



### Docker containers

To get started, install [Docker CE](https://docs.docker.com/install/) or [Docker for Windows](https://docs.docker.com/docker-for-windows/install/).
You can either use already prepared images available at docker hub [![Docker hub](https://img.shields.io/badge/docker-hub-blue.svg?colorA=2271b8&colorB=dc750d&logo=docker&style=flat-square&logoColor=FFF)](https://hub.docker.com/u/flow123d/) for stable versions or use developer images to compile the Flow123d for developer versions.

To start type
```sh
$> docker run flow123d/3.0.1 flow123d --version
00:00:00.035      This is Flow123d, version 3.0.1 commit: 0416d9b
00:00:00.035      Branch: 3.0.1

```

To actually process some file, you probably want to mount a folder, in which is the file located. Assuming the problem file is located in /foo/bar.yaml:
```sh
docker run -ti --rm -euid=$(id -u) -v /foo:/foo flow123d/3.0.1 flow123d --version /foo/bar.yaml
```

Jump to one of these section to learn more:
  - [Linux (stable)](#linux-stable)
  - [Linux (dev)](#linux-dev)
  - [Windows (stable)](#windows-stable)
  - [Windows (dev)](#windows-dev)

#### Linux (stable)

To use specific stable version type the following:
```sh
docker pull flow123d/3.0.1
docker run flow123d/3.0.1 flow123d --version
```
or download the archive from http://flow123d.github.io/, untar the archive and run install.sh located in a bin folder.

#### Linux (dev)

To build the latest version, you can use docker images which
are fully equipped with necessary libraries and frameworks.

Clone the repository, then enter the main directory and copy the docker configuration file:
```sh
git clone https://github.com/flow123d/flow123d.git
cd flow123d/
cp config/config-jenkins-docker-debug.cmake config.cmake
```

After the configuration has been prepared, enter the docker image via utility wrapper
script `fterm`:

> An `fterm` script is a wrapper for the docker. It is recommended to
> use this tool instead of raw docker calls as it automatically mounts necessary paths
> and takes care of various other necessities such as preserving file permission.

```sh
bin/fterm
```

In order to use a release version of the libraries, first, copy a `release` config file:
```sh
cp /config/config-jenkins-docker-release.cmake config.cmake
```

An `fterm` script is by default using a docker image which contains a `debug` version of libraries.
They can be useful to a developer but for production use you may want to switch to a `release` version.
To use a `release` libraries run `fterm` with a flag `rel` (as for a `release`):

```sh
bin/fterm rel
```

This process may take **several minutes** for the first time
(it greatly depends on your internet connection). It will download a debug (or a release)
docker image from the [docker hub](https://hub.docker.com/u/flow123d/) and personalize
the image for the current user.

After the command, you should enter the docker container and see something like this:
```
 ___ _            _ ___ ____    _  
| __| |_____ __ _/ |_  )__ / __| |
| _|| / _ \ V  V / |/ / |_ \/ _  |
|_| |_\___/\_/\_/|_/___|___/\__,_|
                         dbg-3.0.1
jan-hybs@dbg-3.0.1 ~/projects/flow123d/flow123d
```

You should be now all set up for the compilation, simply run `make` while inside
the docker container and wait for the compilation to finish. This may
take up to **15 minutes**, depending on your computer's performance. You can
speed up the compilation by specifying the `-j` flag.


#### Windows (stable)
To use Flow123d under Windows, download your desired version from http://flow123d.github.io/. After that execute installer adn follow the instructions on screen. Installer can download up to 1GB.

#### Windows (dev)

Make sure to have [Docker for Windows](https://docs.docker.com/docker-for-windows/install/) in order to use dev version of Flow123d. The rest should be the same as in section [Linux (dev)](#linux-dev)

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

## Getting started
Please refer to a **User Guide and Input Reference manual** available
at our [official website](http://flow123d.github.io/) where there is a entire section dedicated
to this topic. You can find step-by-step tutorial explaining geometries, `yaml` input files
and more. Below you can see a result from the tutorial problem.
![](/doc/graphics/figure.png)



## Developers
### Troubleshooting

  * When problem occurs during the compilation process it may be due to a leftover files in a build folder.
  Cleaning this directory can solve this issue. You can either remove `build-<branch>` folder
  (the folder is located one level above repository root) via
  `make clean-all`, which removes build folders and also remove any symlinks.  
  To clean all the build folders manually run <code>rm -rf ../build-*</code> while in a repository root.  
  **Running `rm -rf` can quite easily cause a lot of damage, double check that you're
  in a correct folder.**

  * During an installation under Windows, some scenarios can cause problems. Please refer to
  [an installation guide](https://docs.docker.com/toolbox/toolbox_install_windows/) for a
  Docker Toolbox. You can also check out
  [Troubleshooting page](https://docs.docker.com/toolbox/faqs/troubleshoot/) where the most
  common error are described and solved.


### Building the reference manual

The reference manual can be built by while in docker container
```sh
make ref-doc
```
To copy out reference manual from docker use command
[`docker cp`](https://docs.docker.com/engine/reference/commandline/cp/).
