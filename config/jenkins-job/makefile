# Makefile for uploading configuration via jenkins-jobs-builder
# author: Jan Hybs


# echo cmake -DINPUT=main.yaml
# ${CMAKE_COMMAND} -DINPUT=jenkins_jobs_template.yaml -DOUTPUT=${JENKINS_JOBS_FILE} -P ${CMAKE_CURRENT_SOURCE_DIR}/expand_placeholders.cmake

# default configuration for jenkins-jobs-builder
ifdef CONFIG
CONFIG_ARG=--conf $(CONFIG)
endif

# target which replaces placeholders in all *.yaml files in jobs dir
# located in jobs dir
.PHONY: jobs/**/*.yaml jobs/*.yaml
jobs/**/*.yaml jobs/*.yaml:
	cmake -DINPUT=$@ -DOUTPUT=tmp/$@ -P expand_placeholders.cmake


# target which call replacement on all *.yaml files in jobs dir
.PHONY: expand
expand: jobs/**/*.yaml jobs/*.yaml


# target which runs "test" on jenkins-jobs
.PHONY: test
test: expand
	@echo "-----------------------------------------------"
	@echo "Using configuration scripts $(CONFIG)"
	@echo "-----------------------------------------------"
	jenkins-jobs $(CONFIG_ARG) test tmp/jobs/main.yaml


# target which runs "test" on jenkins-jobs
.PHONY: upload
upload: expand
	@echo "-----------------------------------------------"
	@echo "Using configuration scripts $(CONFIG)"
	@echo "-----------------------------------------------"
	jenkins-jobs $(CONFIG_ARG) --ignore-cache update tmp/jobs/main.yaml


# target which uploads configuration via jenkins-jobs to server
.PHONY: delete
delete: #which
	jenkins-jobs $(CONFIG_ARG) --ignore-cache delete-all


#.PHONY: which
#which:
#	@echo "-----------------------------------------------"
#	@echo "Using configuration scripts $(CONFIG)"
#	@echo "-----------------------------------------------"
#	cat $(CONFIG)
#	@echo "-----------------------------------------------"

.PHONY: help
help:
	@echo "Usage:"
	@echo "    - test:   expands placeholders and tests configuration without"
	@echo "              upload and prints result xml"
	@echo ""
	@echo "    - upload: expands placeholders and uploads configuration to server"
	@echo "              by default config file is used."
	@echo "              You can override this behaviour by specifying variable"
	@echo "              CONFIG in make step, example:"
	@echo "              make CONFIG=foo.ini upload"
	@echo ""
	@echo "    - expand: expands placeholders in all files located in jobs folder"
	@echo "              ending with .yaml.template. This action will preserve originals"
	@echo "              and created new files where suffix .template will be removed"
	@echo "              main.yaml.template -> main.yaml"
	@echo ""
	@echo "    - delete: *** deletes all jobs from server ***"
	@echo "              this action cannot be undone"
	@echo ""
	@echo "    - help:   print this message"
