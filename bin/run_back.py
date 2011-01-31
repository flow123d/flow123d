#!/usr/bin/python

# Script for submitting flow123 tasks to PBS in background.
# It accepts path to a directory with esired task data to be run and path to flow_run.sh script.
# It creates a special directory  for each number of processors the task is about to be run on,
# copies the content of the passed directory to all the newly created directories and then
# submits each job to PBS. Then it periodically executes check if all the jobs are finished,
# after all a report is generated.

import sys
import os
import datetime
import shutil
import time
import commands
import re

def main():

    #check if the argument were passed correctly
    try:
    	dir = sys.argv[1]
        flowrun = sys.argv[2]
    except IndexError:
        sys.exit()
        pass

    #create base directory
    now = datetime.datetime.now()
    currentDir = os.getcwd()
    newDir = os.path.join(currentDir, now.strftime("%Y_%m_%d-%H.%M.%S"))

    if os.path.exists(newDir) == 0:
        os.mkdir(newDir)

    measurements = 3
    proc = [1, 2, 4, 6, 8, 10, 12, 16, 20, 24, 32]
    qsubs = {}  #dictionary of created .qsub files and corresponding processor counts
    qsubJobs = []  #list of job-IDs assigned by PBS

    #create directory for each number of processors the task is about to run on
    for nproc in proc:
        targetDir = os.path.join(newDir, str(nproc))
        if os.path.exists(targetDir) == 0:
            os.mkdir(targetDir)

        #copy files from the source directory to the newly created one
        #also try to find either flow.ini or trans.ini file
        iniFile = ""
        fileNames = os.listdir(dir)
        for fileName in fileNames:
            filePath = os.path.join(dir, fileName)

            #copy everything except the .pos files
            if os.path.isfile(filePath) and os.path.splitext(fileName)[1] != ".pos":
                shutil.copy(filePath, targetDir)

            if fileName == "flow.ini" or fileName == "trans.ini":
                iniFile = os.path.join(targetDir, fileName)

        #create .qsub file
        qsub = str(createQsub(newDir, flowrun, iniFile, nproc))
        if len(qsub) > 0:
            qsubs[qsub] = nproc

    for i in range (measurements):
        #submit all the .qsub files
        for qsub in qsubs.iterkeys():
            #for Hydra:
            out = commands.getstatusoutput("qsub -e " + newDir + " -o " + newDir + " -pe orte " + str(qsubs[qsub]) + " " + qsub)
            if out[0] == 0:   #no error
                #for Hydra, we need to filter out the response got in the form 'Your job _number_ (_file_) has been submitted'
                pattern = re.compile("Your job [0-9]+", re.IGNORECASE)
                match = pattern.match(out[1])
                if match:
                    qsubJobs.append(match.group()[9:])

            #for Rex:
##            out = commands.getstatusoutput("qsub -e " + newDir + " -o " + newDir + " " + qsub)
##            if out[0] == 0:   #no error
##                qsubJobs.append(out[1])


    origList = list(qsubJobs)

    while len(qsubJobs) > 0:
        #every 5 minutes, ask the qstat for information about the jobs
        time.sleep(300)

        for jobId in origList:
            out = commands.getstatusoutput("qstat -j " + str(jobId))
            #determine if the job is running/waiting or it has already finished
            if out[0] == 0:
                if len(out[1]) < 100: #FIXME
                    qsubJobs.remove(jobId)

    # all the jobs have been processed
    if len(qsubJobs) == 0:
        os.system("./make_report.py \"" + str(dir) + "\"")



def createQsub(dir, flowrun, iniFile, nproc):
    fileName = os.path.join(dir, "qsub" + str(nproc) + ".qsub")
    try:
        qsubFile = file(fileName, 'w')
        targetDir = os.path.join(dir, str(nproc))

        #for Hydra:
        qsubFile.write("#!/bin/bash\n")
        qsubFile.write("#$ -S /bin/bash\n")
        qsubFile.write("export OMPI_MCA_plm_rsh_disable_qrsh=1\n")
        qsubFile.write("\n")
        qsubFile.write(flowrun + " $NSLOTS " + targetDir + " " + iniFile)

        #for Rex:
##        qsubFile.write("#PBS -S /bin/bash\n")
##        qsubFile.write("#PBS -j oe\n")
##        qsubFile.write("#PBS -N flow_" + str(nproc) + "\n")
##        qsubFile.write("#PBS -m bae\n")
##        qsubFile.write("#PBS -l walltime=10000\n")
##        qsubFile.write("#PBS -l select=1:ncpus=" + str(nproc) + ":host=rex\n")
##        qsubFile.write("#PBS -l place=free:shared\n")
##        qsubFile.write("\n")
##        qsubFile.write(". /opt/intel/Compiler/11.1/046/bin/iccvars.sh ia64\n")
##        qsubFile.write(". /usr/share/modules/init/bash\n")
##        qsubFile.write("module add mpt\n")
##        qsubFile.write("export KMP_MONITOR_STACKSIZE=64K\n")
##        qsubFile.write("\n")
##        qsubFile.write(flowrun " " + str(nproc) + " " targetDir + " " + iniFile)

        qsubFile.close()

        return fileName
    except IOError:
        #do nothing
        pass

    return ""

if __name__ == '__main__':
    main()