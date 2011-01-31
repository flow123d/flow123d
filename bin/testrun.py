#!/usr/bin/python

# Main script for running flow123 tasks.
# It accepts path to a directory with desired task data to be run.
# It starts a background process (run_back.py) which creates a special directory
# for each number of processors the task is about to be run on, copies the content
# of the user-specified directory to all the newly created directories and then
# submits each job to PBS. Then the background process periodically executes
# the qstat command to determine when all the jobs are finished, after all
# a report is generated.

import os
import sys
import subprocess

def main():
    #check whether the user entered an existing directory
    try:
    	dir = sys.argv[1]
    except IndexError:
    	sys.exit("No directory specified")

    dir = os.path.abspath(dir)
    if os.path.exists(dir) == 0 or os.path.isdir(dir) == 0:
    	sys.exit("Specified directory doesn't exist")


    #find the flow_run.sh script
    flowrun = findFlowRun(dir)
    if os.path.exists(flowrun) == 0 or os.path.isfile(flowrun) == 0:
    	sys.exit("Cannot find the flow_run.sh file")

    print "Starting background process..."

    subprocess.Popen([os.path.join(os.getcwd(), "run_back.py"), dir, flowrun])

def findFlowRun(dir):
    flowrun = os.path.join(dir, "flow_run.sh")
    if os.path.exists(flowrun) and os.path.isfile(flowrun):
        return flowrun
    else:
        parentDir = os.path.dirname(dir)
        if os.path.exists(parentDir) and os.path.isdir(parentDir):
            return findFlowRun(parentDir)
        else:
            return ""


if __name__ == '__main__':
    main()
