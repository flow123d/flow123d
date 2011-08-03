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

    print "Starting background process..."

    subprocess.Popen([os.path.join(os.getcwd(), "run_back.py"), dir])


if __name__ == '__main__':
    main()
