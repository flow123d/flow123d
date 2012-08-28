#!/usr/bin/python

import os
import sys
import re

class TimerInfo:
    def __init__(self, proc, taskSize):
        self.taskSize = taskSize
        self.proc = proc
        self.tags = {}
    def setTag(self, tag, time):
        self.tags[tag] = time



#process the output of timers with following tags
tagsToProcess = ["WHOLE PROGRAM","SOLVING MH SYSTEM"]	#,"TRANSPORT"
includeSubdirs = 1

def main():
    #check whether the user entered an existing directory
    try:
    	dir = sys.argv[1]
    except IndexError:
    	sys.exit("No directory specified")

    dir = os.path.abspath(dir)
    if os.path.exists(dir) == 0 or os.path.isdir(dir) == 0:
    	sys.exit("Specified directory doesn't exist")


    #find the longest tag name (for formatting the output)
    longestTag = -1
    for timerTag in tagsToProcess:
        length = len(timerTag)
        if longestTag < length:
            longestTag = length

    #dictionary containing the information gained from the profiler output files,
    #the key of the dictionary is the number of processors, the value is another dictionary
    #containing timer tags as keys and TagInfo objects as values
    fileInfo = []

    processDir(dir, includeSubdirs, fileInfo)

    #check if all the tasks were run on a constant processors count or not
    #and if they differ in the task size or not
    constantProcNumber = 1
    constantTaskSize = 1
    if len(fileInfo) > 0:
        for info in fileInfo:
            if info.proc != fileInfo[0].proc:
                constantProcNumber = 0
            if info.taskSize != fileInfo[0].taskSize:
                constantTaskSize = 0

    if constantProcNumber == 0:
        #sort by the number of processors
        fileInfo.sort(key=lambda x: x.proc)
    elif constantTaskSize == 0:
        #sort by the task size
        fileInfo.sort(key=lambda x: x.taskSize)

    #generate report
    #create the output file
    outFile = file(os.path.join(dir, "report"), "w")

    #the first line
    outFile.write("n. proc".ljust(longestTag + 2))
    for info in fileInfo:
        outFile.write(str(info.proc).rjust(7))

    outFile.write("\n")

    for timerTag in tagsToProcess:
        #write info for each tag we wanted to process
        outFile.write(timerTag.ljust(longestTag + 2))

        if len(fileInfo) > 0:

            #write times on the first line
            for info in fileInfo:
                if info.tags.has_key(timerTag):
                    outFile.write(str(info.tags[timerTag]).rjust(7))
                else:
                    outFile.write("".rjust(7))

            outFile.write("\n")
            outFile.write("".ljust(longestTag + 2))

            if constantProcNumber == 0:
                #compute and write effectivities on the second line
                minProcInfo = fileInfo[0]
                for info in fileInfo:
                    if info.tags.has_key(timerTag) and minProcInfo.tags.has_key(timerTag):
                        effectivity = ((minProcInfo.tags[timerTag] * minProcInfo.proc) / (info.tags[timerTag] * info.proc)) * (info.taskSize / minProcInfo.taskSize)
                        outFile.write(str("%.2f" % effectivity).rjust(7))

                    else:
                        outFile.write("".rjust(7))
            elif constantTaskSize == 0:
                #compute and write effectivities on the second line
                minTaskSizeInfo = fileInfo[0]
                for info in fileInfo:
                    if info.tags.has_key(timerTag) and minProcInfo.tags.has_key(timerTag):

                        effectivity = (minProcInfo.tags[timerTag] / info.tags[timerTag]) * (info.taskSize / minProcInfo.taskSize)
                        outFile.write(str("%.2f" % effectivity).rjust(7))

                    else:
                        outFile.write("".rjust(7))

        outFile.write("\n")

    outFile.close()
    
    graph1(dir,timerTag,longestTag,fileInfo)
    graph2(dir,timerTag,longestTag,fileInfo,constantProcNumber,constantTaskSize)
    
    
    sys.exit()
    
    
    ###
    
    # T(p)
    
def graph1(dir,timerTag,longestTag,fileInfo):
    
    for timerTag in tagsToProcess:         
        outFile = file(os.path.join(dir, 'graph.dat'), "w")
        outFile.write("# n. proc".ljust(2))
	for i in range(len(tagsToProcess)):  
         outFile.write(tagsToProcess[i].rjust(longestTag+2))
        outFile.write("\n")
        for info in fileInfo:
         outFile.write(str(info.proc).ljust(2))
         #maxtime = 0.0
         for timerTag in tagsToProcess: 
          outFile.write(str(info.tags[timerTag]).rjust(longestTag+3))
          #if info.tags[timerTag] > maxtime:
          # maxtime = info.tags[timerTag]
         outFile.write("\n")  
        outFile.close()
    for i in range(len(tagsToProcess)):    
     outFile = file(os.path.join(dir, tagsToProcess[i]+".gps",), "w")
     outFile.write("#!/usr/bin/gnuplot\n")
      #outFile.write('set title "'+tagsToProcess[i]+'"\n')
     outFile.write('set yrange [0:'+str()+']\n')	
     #outFile.write('set xrange [1:2]\n') len(fileInfo) fileInfo.__len__()-1
     outFile.write('set xrange ['+str(fileInfo[0].proc).rjust(0)+':'+str(fileInfo[len(fileInfo)-1].proc).rjust(0)+']\n')
     outFile.write('set xtics 1, 1\n')
     outFile.write('set mxtics 1\n') #linespoints boxes
     outFile.write("plot '"+os.path.join(dir, "graph.dat")+"' using 1:"+str(i+2)+" title '"+tagsToProcess[i]+"' with linespoints\n")
     outFile.write("pause -1")
     outFile.close()
     os.system("gnuplot '"+os.path.join(dir, tagsToProcess[i])+".gps'")
     
     
     ###
     
     # E(p), size = const
     
     
def graph2(dir,timerTag,longestTag,fileInfo,constantProcNumber,constantTaskSize):     

    for timerTag in tagsToProcess:         
        outFile = file(os.path.join(dir, 'graphx.dat'), "w")
        outFile.write("# n. proc".ljust(2))
	for i in range(len(tagsToProcess)):  
         outFile.write(tagsToProcess[i].rjust(longestTag+2))
        outFile.write("\n")
        
     #   for info in fileInfo:
     #    outFile.write(str(info.proc).ljust(2))
       #  #maxtime = 0.0
       #  for timerTag in tagsToProcess: 
       #   outFile.write(str(info.tags[timerTag]).rjust(longestTag+3))
          #if info.tags[timerTag] > maxtime:
          # maxtime = info.tags[timerTag]
       #  outFile.write("\n")  
         
    #for timerTag in tagsToProcess:
    #for info in fileInfo:  
    if constantProcNumber == 0:
                #compute and write effectivities on the second line
                minProcInfo = fileInfo[0]
                for info in fileInfo:
                    outFile.write(str(info.proc).ljust(2))
                    for timerTag in tagsToProcess:
                     if info.tags.has_key(timerTag) and minProcInfo.tags.has_key(timerTag):
                        effectivity = ((minProcInfo.tags[timerTag] * minProcInfo.proc) / (info.tags[timerTag] * info.proc)) * (info.taskSize / minProcInfo.taskSize)
                        outFile.write(str("%.2f" % effectivity).rjust(7))
                        #outFile.write("\n")
                     else:
                        outFile.write("".rjust(7))
                    outFile.write("\n")   
    elif constantTaskSize == 0:
                #compute and write effectivities on the second line
                minTaskSizeInfo = fileInfo[0]
                #for info in fileInfo:
                for info in fileInfo:  
                    outFile.write(str(info.proc).ljust(2))
                    for timerTag in tagsToProcess:
                     if info.tags.has_key(timerTag) and minProcInfo.tags.has_key(timerTag):
                        effectivity = (minProcInfo.tags[timerTag] / info.tags[timerTag]) * (info.taskSize / minProcInfo.taskSize)
                        outFile.write(str("%.2f" % effectivity).rjust(7))   
                     else:
                        outFile.write("".rjust(7))
                    outFile.write("\n")   
     
    outFile.close()
    for i in range(len(tagsToProcess)):    
     outFile = file(os.path.join(dir, tagsToProcess[i]+".gpse",), "w")
     outFile.write("#!/usr/bin/gnuplot\n")
      #outFile.write('set title "'+tagsToProcess[i]+'"\n')
     outFile.write('set yrange [0:'+str()+']\n')	
     #outFile.write('set xrange [1:2]\n') len(fileInfo) fileInfo.__len__()-1
     outFile.write('set xrange ['+str(fileInfo[0].proc).rjust(0)+':'+str(fileInfo[len(fileInfo)-1].proc).rjust(0)+']\n')
     outFile.write('set xtics 1, 1\n')
     outFile.write('set mxtics 1\n') #linespoints boxes
     outFile.write("plot '"+os.path.join(dir, "graphx.dat")+"' using 1:"+str(i+2)+" title '"+tagsToProcess[i]+"' with linespoints\n")
     outFile.write("pause -1")
     outFile.close()
     os.system("gnuplot '"+os.path.join(dir, tagsToProcess[i])+".gpse'")     
       
       
     ###  
     
     # E(p), size != const
      
    
    
      #Create a title:                  > set title "Force-Deflection Data" 
      #Put a label on the x-axis:       > set xlabel "Deflection (meters)"
      #Put a label on the y-axis:       > set ylabel "Force (kN)"
      #Change the x-axis range:         > set xrange [0.001:0.005]
      #Change the y-axis range:         > set yrange [20:500]
    
    

def processDir(dir, includeSubdirs, fileInfo):

    patternNumproc = re.compile("No\. of processors: (?P<num>[0-9]+)")
    patternTaskSize = re.compile("Task size: (?P<num>[0-9]+)")


    fileNames = os.listdir(dir)
    for fileName in fileNames:
        filePath = os.path.join(dir, fileName)

        if os.path.isdir(filePath) and includeSubdirs:
            processDir(filePath, includeSubdirs, fileInfo)
        elif os.path.isfile(filePath) and filePath.endswith("out"):
            #read content of the file
            f = file(filePath, "r")
            fileContent = f.read()
            f.close()

            #find the number of processors and task size using the regular expressions
            numprocmatch = patternNumproc.search(fileContent)
            tasksizematch = patternTaskSize.search(fileContent)

            if numprocmatch and tasksizematch:
                size = int(tasksizematch.group("num"))
                proc = int(numprocmatch.group("num"))
                for timerTag in tagsToProcess:
                    #for each tag we want to process, create the regular expression
                    patternTag = re.compile("\s*"+timerTag+"\s+[0-9]+\s+(?P<num>[0-9]+(\.[0-9]+)?)")
                    tagmatch = patternTag.search(fileContent)

                    if tagmatch:
                        time = float(tagmatch.group("num"))

                        found = 0
                        for info in fileInfo:
                            if info.proc == proc and info.taskSize == size:
                                found = 1
                                if info.tags.has_key(timerTag):
                                    if info.tags[timerTag] > time:
                                        info.tags[timerTag] = time
                                else:
                                    info.setTag(timerTag, time)
                                break
                        if found == 0:
                            info = TimerInfo(proc, size)
                            info.setTag(timerTag, time)
                            fileInfo.append(info)

if __name__ == '__main__':
    main()
