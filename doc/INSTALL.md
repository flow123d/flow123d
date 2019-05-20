# Installation

This guide summarises how to install Flow123d. They are two ways
to install depending on whether you want to use **stable**
or **developer** version of Flow123d:

  1) [**Stable version**](#stable-version) using prebuilt Docker images:

     This is the easiest and the most reliable way to start using Flow123d.
     It is available for Linux and Windows users.
     Windows platform support is limited and requires
     [Docker for Windows](https://docs.docker.com/docker-for-windows/install/).
     Linux platform requires [Docker CE](https://docs.docker.com/install/).

  2) [**Dev version**](#dev-version) i.e. building Flow123d from source:

      This option requires a moderate to high understanding
      of a operating systems but offers mostly higher performance.

      1) with usage of container technology
      2) without any containerization, i.e. native Linux

**For installation in high performance computing clusters we recommend always building from source**.

## Stable version

### Linux (stable)
All prebuilt Docker images are hosted at [![Docker hub](https://img.shields.io/badge/docker-hub-blue.svg?colorA=2271b8&colorB=dc750d&logo=docker&style=flat-square&logoColor=FFF)](https://hub.docker.com/u/flow123d/).
You can `docker pull` any of them, but to avoid mounting issues,
you can use following bash snippet, which will download `flow123d`
script allowing you to explore software further.

In a terminal simply run:
```sh
curl -s https://flow.nti.tul.cz/get | bash
flow123d
```

And you should see welcome screen. For the first time,
it may take a minute to download Docker Image (your internet
connection will do all the heavy lifting).
```
 ___ _            _ ___ ____    _  
| __| |_____ __ _/ |_  )__ / __| |
| _|| / _ \ V  V / |/ / |_ \/ _  |
|_| |_\___/\_/\_/|_/___|___/\__,_|
                        3.0.1
flow@flow:3.0.1 /home/jan-hybs
```

---

You can also download archive from [official site](http://flow123d.github.io/) which when unzipped contains
`install.sh` which, when executed, will download image and
allows you to do the same.

```sh
wget http://flow.nti.tul.cz/packages/2.2.1_release/flow123d_2.2.1_linux_install.tar.gz
tar -xvf flow123d_2.2.1_linux_install.tar.gz
flow123d_2.2.1/install.sh
```

### Windows (stable)

To use Flow123d under Windows, download your desired version from [official site](http://flow123d.github.io/). After that execute installer and follow the instructions on screen. Installer can download up to 1GB.


## From source (dev version)


### Linux (dev)
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


### Windows (dev)

Make sure to have [Docker for Windows](https://docs.docker.com/docker-for-windows/install/) in order to use dev version of Flow123d. The rest should be the same as in section [Linux (dev)](#linux-dev)

## From source (native)
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
[flow123d-docker-images](https://github.com/flow123d/flow123d-docker-images)
repository, where you can learn more.

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
