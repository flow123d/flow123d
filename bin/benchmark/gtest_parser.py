import sys
import os
import xml.etree.ElementTree as ET
import platform
import re
import tempfile
import subprocess

class ExcWrongTag(Exception):
    def __init__(self, value):
        self.value=value
    def __str__(self):
        return "Wrong tag: " + repr(self.value)
      
class ExcMessage(Exception):
    def __init__(self, value):
        self.value=value
    def __str__(self):
        return repr(self.value)      
        
def make_node(type):
    return { 'attrib' : { 'type' : type }, 'childs' : {} }

def parse_gtest_xml( filename ):
    """
    Returns dictionary with data from a google test JUnit XML report.
    Dictionary structure:
    { attrib: { type: "test" }
      childs: {
          suite1 : {
            attrib: { type: "suite", disabled: false }
            childs: { 
              case1: { 
                attrib: {type: "case", run: true, time: 123  }
              }}}}}
     
    """
    tree = ET.parse(filename)
    testsuits = tree.getroot()
    
    if (testsuits.tag != "testsuites"):
        raise ExcWrongTag(testsuits.tag)
   
    result=make_node("test")
    for suite in testsuits:
        if (suite.tag != "testsuite"):
            raise ExcWrongTag(suite.tag)
        
        r_suite=result['childs'][suite.attrib['name']]=make_node("suite")        
        r_suite['attrib']['disabled']= ( suite.attrib['disabled'] != 0 )
        for case in suite:
            if (case.tag != "testcase"):
                raise ExcWrongTag(case.tag)
            
            r_case=r_suite['childs'][case.attrib['name']]=make_node("case")
            r_case['attrib']['time']=case.attrib['time']
            r_case['attrib']['run']= ( case.attrib['status'] == 'run' )
    return result        

def get_host_info():
    """
    Returns a dictionary with information about testrun environment, namely
    specification of host and code version
    {
      os:               # {win|linux}_{32|64}
      hostname:         # local host name, address 
      uname:
      frequency:        # in MHz
      benchmark:        # time of sample code snipet
    """
    result={}
    uname=platform.uname()
    result['uname']={
      'system':uname[0],
      'hostname':uname[1],
      'release':uname[2],
      'version':uname[3],
      'machine':uname[4],
      'processor':uname[5]}
    if (uname[0] == "Linux") :
        subprocess.check_output(['cat','/proc/cpuinfo'])

    # frequency from system  
    result['frequency']=None
    if (uname[0] == "Windows") :
        import wmi
        
        c=wmi.WMI()
        result['frequency']=float(c.Win32_Processor()[0].MaxClockSpeed)
    else :
        out = subprocess.check_output(['cat','/proc/cpuinfo'])
        match=re.search('cpuMHz\s*:      ([0-9.]*)',out)
        if (match) :
            result['frequency']=float(match.group(1))
    
    # simple speed test
    import timeit
    time=timeit.Timer('for i in range(0,100) : x+=x*list[x%1000];','x=1;list=range(1,1000)').timeit(100)
    print time
    result['benchmark']=time
        
    result['os']=uname[0]+"_"+uname[4]
    result['hostname']=platform.node() 
    
    return result


def get_source_info():
    """
      git_branch:       
      git_remote:       # URL of remote of actual branch
      git_commit:       # hash
      git_describe:     # 
      git_message:      # commit message first line
    """
    result={}
    
    out=subprocess.check_output(['git', 'branch', '--list', '-vv'])
    match=re.search('\* (\w*)\s*([0-9a-f]*) (\[(.*)\] )?(.*)\n',out)    
    result['git_branch']=match.group(1)
    result['git_short_commit']=match.group(2)
    remote_branch=match.group(4)
    result['git_message']=match.group(5)
    
    print remote_branch
    if (remote_branch != None) :
        remote=re.match('(\w*)/(.*)', remote_branch).group(1)
        out=subprocess.check_output(['git','remote','-v'])
        result['git_remote']=re.search('('+remote+')\s*([^ ]*) \(fetch\)\n', out).group(2)
    else :
        result['git_remote']=None
    
    result['git_describe']=subprocess.check_output(['git','describe'])  
    return result
  
def get_test_data():
    pass
  
def send_to_server(report):
    import json
    import requests
    
    url='http://localhost:8000'
    headers={'content-type': 'application/json'}
    r=requests.post(url, data=json.dumps(report),headers=headers)
    if (r.status_code != 201) :
        raise(ExcMessage("Server do not accept the report."))
  
def make_report():
    report={
      'host_info' : get_host_info(),
      'source_info' : get_source_info(),
      'test_data' : get_test_data()
    }  
      
    return report
# Main -----------------------------------

send_to_server(make_report())

