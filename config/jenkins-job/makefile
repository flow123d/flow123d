# Makefile for uploading configuration via jenkins-jobs-builder 
# author: Jan Hybs


# echo cmake -DINPUT=main.yaml
# ${CMAKE_COMMAND} -DINPUT=jenkins_jobs_template.yaml -DOUTPUT=${JENKINS_JOBS_FILE} -P ${CMAKE_CURRENT_SOURCE_DIR}/expand_placeholders.cmake

# default configuration for jenkins-jobs-builder
CONFIG=jenkins.ini

# target which replaces placeholders in all *.yaml.template files
# located in jobs dir
.PHONY: jobs/*.yaml.template
jobs/*.yaml.template:
	cmake -DINPUT=$@ -DOUTPUT=$(@:%.template=%) -P expand_placeholders.cmake


# target which call replacement on all yaml.template files
.PHONY: expand
expand: jobs/*.yaml.template


# target which runs "test" on jenkins-jobs
.PHONY: test
test: expand
	@echo "-----------------------------------------------"
	@echo "Using configuration scripts $(CONFIG)"
	@echo "-----------------------------------------------"
	jenkins-jobs --conf $(CONFIG) test jobs/main.yaml


# target which runs "test" on jenkins-jobs
.PHONY: upload
upload: expand
	@echo "-----------------------------------------------"
	@echo "Using configuration scripts $(CONFIG)"
	@echo "-----------------------------------------------"
	jenkins-jobs --conf $(CONFIG) --ignore-cache update jobs/main.yaml


# target which uploads configuration via jenkins-jobs to server
.PHONY: delete
delete: which
	jenkins-jobs --conf $(CONFIG) --ignore-cache delete-all


.PHONY: which
which:
	@echo "-----------------------------------------------"
	@echo "Using configuration scripts $(CONFIG)"
	@echo "-----------------------------------------------"
	cat $(CONFIG)
	@echo "-----------------------------------------------"

.PHONY: help
help:
	@echo "Usage:"
	@echo "    - test:   expands placeholders and tests configuration without"
	@echo "              upload and prints result xml"
	@echo ""
	@echo "    - upload: expands placeholders and uploads configuration to server"
	@echo "              by default $(CONFIG) is used."
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
	