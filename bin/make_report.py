#!/usr/bin/env python

import os
import sys
import re

def main():
    #check whether the user entered an existing directory
    try:
    	dir = sys.argv[1]
    except IndexError:
    	sys.exit("No directory specified")

    dir = os.path.abspath(dir)
    if os.path.exists(dir) == 0 or os.path.isdir(dir) == 0:
    	sys.exit("Specified directory doesn't exist")

    patternPath = re.compile("Path: .+")
    patternNumproc = re.compile("n\. proc: (?P<num>[0-9]+)")
    patternNumits = re.compile("Lin Solver: its: (?P<num>[0-9]+)")
    patternSolving = re.compile("Last. SOLVING MH SYSTEM; period: [0-9]+\.[0-9]+ sec; total:  (?P<total>[0-9]+\.[0-9]+) sec\.")
    patternTotal = re.compile("Last. WHOLE PROGRAM; period: [0-9]+\.[0-9]+ sec; total:  (?P<total>[0-9]+\.[0-9]+) sec\.")
    patternMatmult = re.compile("^MatMult.+ (?P<num>[0-9]+)", re.MULTILINE)
    patternKSPSolve = re.compile("^KSPSolve.+ (?P<num>[0-9]+)", re.MULTILINE)

    #dictionary containing number of processors as the key and information gained from the output file as a tuple in value
    #value has the form (numits, solvingtime, totaltime, matMultMFLOPS, KSPSolveMFLOPS)
    fileInfo = {}
    fileNames = os.listdir(dir)
    path = ""
    minProcCount = -1
    for fileName in fileNames:
        filePath = os.path.join(dir, fileName)

        if os.path.isfile(filePath):
            #read content of the file
            f = file(filePath, "r")
            fileContent = f.read()
            f.close()

            #find the info about number of processors and solving time
            solvingmatch = patternSolving.search(fileContent)
            totalmatch = patternTotal.search(fileContent)
            numprocmatch = patternNumproc.search(fileContent)
            numitsmatch = patternNumits.search(fileContent)
            matmultmatch = patternMatmult.search(fileContent)
            kspsolvematch = patternKSPSolve.search(fileContent)
            if len(path) == 0:
                m = patternPath.search(fileContent)
                if m:
                    path = m.group()

            if solvingmatch and totalmatch and numprocmatch and numitsmatch and matmultmatch and kspsolvematch:
                solvingtime = float(solvingmatch.group("total"))   #get time
                proc = int(numprocmatch.group("num"))
                if minProcCount < 0:
                    minProcCount = proc
                else:
                    minProcCount = min(minProcCount, proc)

                found = 0
                if fileInfo.has_key(proc):
                    #any file for this number of processors has been olready found. Choose the file with lower SolvingTime
                    if fileInfo[proc][1] < solvingtime:
                        found = 1
                if found == 0:
                    numits = int(numitsmatch.group("num"))
                    totaltime = float(totalmatch.group("total"))
                    matMult = int(matmultmatch.group("num"))
                    kspSolve = int(kspsolvematch.group("num"))
                    fileInfo[proc] = (numits, solvingtime, totaltime, matMult, kspSolve)
                    pass

    #compute effectivities
    # (effectivity for 1 iteration, total eff, matMult eff, KSPsolve eff)
    effectivity = {}
    if minProcCount > 0:
        minProcInfo = fileInfo[minProcCount]
        for proc in fileInfo.keys():
            info = fileInfo[proc]
            eff1iter = float(minProcInfo[0]*minProcInfo[1])/(info[0]*info[1])
            totalEff = float(minProcInfo[1])/info[1]/proc
            matMultEff = (float(minProcInfo[0])/minProcInfo[3])*(info[3]/info[0])
            kspSolveEff = float(info[4])/(proc*minProcInfo[4])
            effectivity[proc] = (eff1iter, totalEff, matMultEff, kspSolveEff)


    #generate report
    #probably slightly unoptimal, but we want to write the values into columns, not into rows
    outFile = file(os.path.join(dir, "result"), "w")
    outFile.write("Directory " + path + "\n")
    outFile.write("---------------------------------------------------\n")
    outFile.write("n. proc".ljust(16))
    for proc in sorted(fileInfo.keys()):
        outFile.write(str(proc).rjust(7))

    outFile.write("\n")

    captions = ("numits", "solvingtime", "totaltime", "matMultMFLOPS", "KSPSolveMFLOPS", "1 iter.", "total", "MatMult", "KSPsolve")
    for i in range(5):
        outFile.write(captions[i].ljust(16))
        for proc in sorted(fileInfo.keys()):
            if type(fileInfo[proc][i]) == float:
                outFile.write(str("%.2f" % fileInfo[proc][i]).rjust(7))
            else:
                outFile.write(str(fileInfo[proc][i]).rjust(7))
        outFile.write("\n")

    outFile.write("----------------effectivity------------------------\n")

    for i in range(5, 9):
        outFile.write(captions[i].ljust(16))
        for proc in sorted(effectivity.keys()):
            outFile.write(str("%.3f" % effectivity[proc][i-5]).rjust(7))
        outFile.write("\n")

    outFile.close()

if __name__ == '__main__':
    main()
