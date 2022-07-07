# This makefile will create docker images from a running container containing
# Flow123d (assuming the project is already compiled)
#
# Variables passed in:
#   container_name=cont{environment}{build-type} 		 # the source container where the Flow123d has been build.
#   install_image="install-{environment}:$imagesversion" 	 # the base image for the installation image
#   target_image=flow123d-{environment}                          # name of the installation image, ci images use: ci-{enviroment}
#   image_tag={flow123d_version}				 # tag of the image, the release tag for the release images, ci images use: master-{commit} 
#   destination=$destination 					 # working directory to prepare the installation 
#

container_name?=contgnurelease
install_image?=install-gnu

target_image?=ci-gnu
destination=$(PWD)/publish_gnu
flow_repo_location=/opt/flow123d/flow123d

# ?= setting of a variable only if it doesn't have a value.

# These variables should not be passed in but need to check Jenkins.
flow_install_location   ?= /opt/flow123d
flow_repo_location      ?= /opt/flow123d/flow123d
uid 				 ?= $(shell id -u)
gid 				 ?= $(shell id -u)
# The 'docker exec' command runs a new command in a running container
flow_version ?= $(strip $(shell docker exec $(container_name) cat $(flow_repo_location)/version ))
git_hash     ?= $(strip $(shell docker exec $(container_name) sh -c "cd $(flow_repo_location) && git rev-parse --short HEAD"))

# default image tag formed from the git branch name and the commit
image_tag   ?= $(strip $(shell docker exec $(container_name) sh -c "cd $(flow_repo_location) && git rev-parse --abbrev-ref HEAD"))_$(git_hash)


# name of the archives from CMake CPack tool
cmake_package_name=Flow123d-$(flow_version)-Linux.tar.gz

# final archive names
base_name=flow123d_$(flow_version)
docker_arch_name=$(base_name)_docker_image.tar.gz
docker_geomop_arch_name=$(base_name)_docker_geomop_image.tar.gz
lin_arch_name=$(base_name)_linux_install.tar.gz
win_arch_name=$(base_name)_windows_install.zip
win_geomop_arch_name=$(base_name)_windows_geomop_install.zip

# current date in two forms
current_date=$(shell date +"%d-%m-%Y %T")
build_date=$(shell date -u +"%Y-%m-%dT%H:%M:%SZ")

# path to the generated pdf
pdf_location=doc/reference_manual/flow123d_doc.pdf
ist_location=doc/reference_manual/input_reference.json

# Commands assignment
dexec=docker exec $(container_name) # execute command which will follow
dcp=docker cp $(container_name)   # Copy files/folders between a container and the local filesystem

target_tagged=flow123d/$(target_image):$(image_tag)
#geomop_image=flow123d/geomop:$(image_tag)_$(git_hash)

help:
	@echo "usage: "
	@echo "  all           packs entire flow"
	@echo "  help          prints this message"
	@echo "  info          prints package configuration"
	@echo "  push-to-hub   pushes docker image $(target_tagged) to the docker hub"


#$(destination)/$(win_arch_name) $(destination)/$(lin_arch_name)
# 'all' is the default target. As in running just 'make' should do the same as 'make all'
all: info linux windows
	ls -la $(destination)


# human readable shortcut
linux: $(destination)/$(lin_arch_name)

# human readable shortcut
windows: $(destination)/$(win_arch_name)


# print package info
info:
	@echo "current configuration"
	@echo "  install_image:$(install_image)"
	@echo "  flow_version: $(flow_version)"
	@echo "  destination:  $(destination)"
	@echo "  git_hash:     $(git_hash)"
	@echo "  image_tag:   $(image_tag)"


# removes unwanted files before publishing to the flow.nti.tul.cz
remove-unwanted:
	rm -rf $(destination)/$(docker_arch_name)
	rm -rf $(destination)/tests
	rm -rf $(destination)/doc
	rm -rf $(destination)/config
	rm -rf $(destination)/bin
	rm -rf $(destination)/win
	rm -rf $(destination)/nsis
	rm -rf $(destination)/install.nsi
	@echo "Following files will be included in package: "
	ls -la $(destination)

# create installation package for windows users
$(destination)/$(win_arch_name): $(destination)/$(docker_arch_name)
	cp -r project/src/windows/* $(destination)/
	echo "$(flow_version)" > $(destination)/version
	echo "$(target_tagged)" > $(destination)/imagename

	# The 'docker run' command first creates a writeable container layer over the specified image, and then starts it using the specified command
	docker run -i --rm -u $(uid):$(gid) -v $(destination):/nsis-project hp41/nsis /nsis-project/install.nsi
	echo "{\"build\": \"$(current_date)\", \"hash\": \"$(git_hash)\"}" > $(destination)/flow123d_$(flow_version)_windows_install.json


# create TAR.GZ package for linux users
$(destination)/$(lin_arch_name): $(destination)/$(docker_arch_name)

	mkdir -p install-linux
	#  mkdir -p make parent directories as needed
	cd install-linux && cmake \
		-DFLOW_VERSION="$(flow_version)" \
		-DFLOW123D_ROOT="$(flow_repo_location)" \
		-DIMAGE_TAG="$(target_tagged)" \
		-DIMAGE_NAME="$(docker_arch_name)" \
		-DDEST="$(destination)" \
		../project

	make -C install-linux package
	mv install-linux/$(base_name).tar.gz $(destination)/$(lin_arch_name)
	echo "{\"build\": \"$(current_date)\", \"hash\": \"$(git_hash)\"}" > $(destination)/$(lin_arch_name:.tar.gz=.json)

ci_docker_file_path=./dockerfile
# create release docker image with flow123d binary and depdendencies
# from docker image flow123d/install
$(destination)/$(docker_arch_name): $(destination)/$(cmake_package_name)
	-@docker rmi -f $(target_tagged) 
	docker commit $(container_name) flow123d/temporary_build
	docker pull flow123d/$(install_image)
	#cp $(destination)/$(cmake_package_name) project/src/docker/create/default/$(cmake_package_name)
	# 'docker build' build an image from a Dockerfile
	docker build \
	 --build-arg base_image=flow123d/$(install_image) \
         --build-arg source_image=flow123d/temporary_build \
         --build-arg source_location=$(flow_install_location)/build_package \
         --build-arg flow_version=$(flow_version) \
         --build-arg flow_install_location=$(flow_install_location) \
         --build-arg git_hash=$(git_hash) \
         --build-arg build_date=$(build_date) \
         --tag $(target_tagged) \
         $(ci_docker_file_path)
	# 'docker save' Save one or more images to a tar archive (streamed to STDOUT by default)
	# docker save $(target_tagged) > $(destination)/$(docker_arch_name)
	#rm -rf project/src/docker/create/default/$(cmake_package_name)

# geomop_docker_file_path=../dockerfiles/geomop-X
# geomop_image: $(destination)/$(docker_arch_name)
# 	docker build \
# 	 --build-arg base_image=$(target_tagged) \
#          --build-arg flow_version=$(flow_version) \
#          --build-arg flow_install_location=$(flow_install_location) \
#          --build-arg git_hash=$(git_hash) \
#          --build-arg build_date=$(build_date) \
#          --tag ${geomop_image} \
#          $(geomop_docker_file_path)
	
# create package from flow123d project
$(destination)/$(cmake_package_name): $(destination)/tests $(destination)/flow123d_$(flow_version)_doc.pdf
	$(dexec)  /bin/bash -c ". $(flow_repo_location)/bin/setvars.sh && make -C $(flow_repo_location) package"
	$(dexec) cp -r $(flow_repo_location)/build_tree/_CPack_Packages/Linux/TGZ/Flow123d-$(flow_version)-Linux $(flow_install_location)/build_package
	# we only use temporary installation and copy it directly form the build image to the taget image
	# TODO: need to copy the package out of the vmount, but only for testing
	
# copy tests outside for system dependent packages, tests and doc do not go into the docker image
$(destination)/tests:
	# create destination folders
	mkdir -p $(destination)/tests
	# clean and copy tests folder

	$(dexec) /bin/bash -c ". $(flow_repo_location)/bin/setvars.sh && make -C $(flow_repo_location) clean-tests"
	docker cp $(container_name):$(flow_repo_location)/tests/. $(destination)/tests

	# delete runtest because we will have to create other runtest for docker
	rm -rf $(destination)/tests/runtest


# call ref doc if pdf was not already created
$(destination)/flow123d_$(flow_version)_doc.pdf:
	$(dexec) /bin/bash -c ". $(flow_repo_location)/bin/setvars.sh && make -C $(flow_repo_location) all"                         # compile just in case
	$(dexec) /bin/bash -c ". $(flow_repo_location)/bin/setvars.sh && make -C $(flow_repo_location) FORCE_DOC_UPDATE=1 ref-doc"   # generate latex doc
	$(dexec) /bin/bash -c ". $(flow_repo_location)/bin/setvars.sh && make -C $(flow_repo_location) html-doc"                    # generate html doc
	$(dexec) /bin/bash -c ". $(flow_repo_location)/bin/setvars.sh && make -C $(flow_repo_location) doxy-doc"                     # generate source doc

	mkdir -p $(destination)/htmldoc
	mkdir -p $(destination)/doxygen
	mkdir -p $(destination)/config/docker/
	mkdir -p $(destination)/bin/

	docker cp $(container_name):$(flow_repo_location)/build_tree/$(pdf_location)             $(destination)/flow123d_$(flow_version)_doc.pdf
	docker cp $(container_name):$(flow_repo_location)/build_tree/htmldoc/html/src/.          $(destination)/htmldoc
	docker cp $(container_name):$(flow_repo_location)/build_tree/doc/online-doc/flow123d/.   $(destination)/doxygen
	docker cp $(container_name):$(flow_repo_location)/$(ist_location)                        $(destination)/input_reference.json
	docker cp $(container_name):$(flow_repo_location)/bin/fterm                              $(destination)/bin/fterm


# will push the images to the hub
# you must be logged in already in order to push the images to docker hub
push-to-hub: $(destination)/$(docker_arch_name)
	-docker push $(target_tagged)
	#-docker push $(geomop_image)
