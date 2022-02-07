# Flow123d &middot; [![Jenkins release](http://ciflow.nti.tul.cz/status/release/Flow123d-ci2runner-release-multijob/4)](http://ciflow.nti.tul.cz:8080/job/Flow123d-gnu-release-multijob/) [![Jenkins debug](http://ciflow.nti.tul.cz:8080/status/debug/Flow123d-gnu-debug-multijob/4)](http://ciflow.nti.tul.cz:8080/job/Flow123d-gnu-debug-multijob/) [![Coveralls master](https://img.shields.io/coveralls/github/flow123d/flow123d.svg?style=flat-square&label=coverage)](https://coveralls.io/github/flow123d/flow123d) [![Docker hub](https://img.shields.io/badge/docker-hub-blue.svg?colorA=2271b8&colorB=dc750d&logo=docker&style=flat-square&logoColor=FFF)](https://hub.docker.com/u/flow123d/) [![CI-HPC](https://img.shields.io/badge/ci--hpc-performace-green.svg?style=flat-square)](http://hybs.nti.tul.cz/ci-hpc/)

*Processes in Fractured Media*

Flow123d is a simulator of underground water flow, transport, and mechanical processes in fractured
porous media. Novelty of this software is support of computations on complex
meshes consisting of simplicial elements of different dimensions. Therefore
we can combine continuum models and discrete fracture network models.
For more information see the project pages:
[flow123d.github.io](http://flow123d.github.io). 


## Getting started
Please refer to a **User Guide and Input Reference manual** available
at our [official website](http://flow123d.github.io/) where there is a entire section dedicated
to this topic. You can find step-by-step tutorial explaining geometries, `yaml` input files
and more. Below you can see a result from the tutorial problem.
![](/doc/graphics/figure.png)


## Installation
You have several options when it comes to installation. You can use
our prebuilt Docker images or build Flow123d from the source (requires a moderate experience with Linux OS).

For detailed instructions, see the [installation guide](doc/INSTALL.md).


## Developers

## Build
Two step build:
```
     host> bin/fterm     # start the docker developing container
container> make all      # produce "build_tree/bin/flow123d"
container> bin/flow123d  # lunch the simulator using a simple wrapper script
```

That would start the debug docker image and build the debug version of the simulator.
For the release (optimized) version run:
```
bin/fterm rel
make all
```

For details see [installation guide](doc/INSTALL.md) or the manual.


### Troubleshooting

  * When problem occurs during the compilation process it may be due to a leftover files in a build folder.
  Cleaning this directory can solve this issue. You can either remove `build-<branch>` folder
  (the folder is located one level above repository root) via
  `make clean-all`, which removes build folders and also remove any symlinks.  
  To clean all the build folders manually run <code>rm -rf ../build-*</code> while in a repository root.  
  **Running `rm -rf` can quite easily cause a lot of damage, double check that you're
  in a correct folder.**

  * The build tools may fail if the root path contains folders with spaces.
  
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


## Singularity

[Singularity](https://singularity.hpcng.org/) is a container system targeting HPC applications. Singularity containers can be created from the docker images and in contrast to the docker one is not allowed to modify system of running containers. 

For large images may be necessary set temp directory.
```sh
export SINGULARITY_TMPDIR="/some_absolute_path/tmp"
```

### Run simulator in singularity container
```sh
singularity exec docker://flow123d/3.1.0 flow123d simulation.yaml
```

### Parallel run, mpiexec out of the image
```sh
module add mpich-3.0.2-gcc
mpiexec -host host1,host2 -np 4 singularity exec docker://flow123d/3.1.0 flow123d simulation.yaml
```

### Build
```sh
     host> git clone https://github.com/flow123d/flow123d.git # clone flow123d repository
     host> singularity shell -B flow123d/:/flow123d docker://flow123d/flow-dev-gnu-rel:3.1.0 # starts developing container
container> cd /flow123d
container> make all # produce "build_tree/bin/flow123d"
container> bin/flow123d  # lunch the simulator
```

### Docker images
Hosted on `dockerhub` under organisation [`flow123d`](https://hub.docker.com/orgs/flow123d/repositories).

**production images** 

tagged by the release version (e.g. 3.1.0)

- [`flow123-gnu`](https://hub.docker.com/repository/docker/flow123d/flow123d-gnu) 
  : based on gnu libraries and tools
- `flow123-intel` : based on intel libraries (intelmpi, mkl, intel compilers)

**alfa images**

built by CI, tagged by `{branch}-{commit}`

- [`ci-gnu`](https://hub.docker.com/repository/docker/flow123d/ci-gnu)
- `ci-intel`

