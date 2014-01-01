import BaseHTTPServer
import pymongo
import json
import dateutil.parser

class ExcAssert(Exception):
    def __init__(self, value):
        self.value=value
    def __str__(self):
        return repr(self.value) 


class ExcMessage(Exception):
    def __init__(self, value):
        self.value=value
    def __str__(self):
        return repr(self.value) 
      
class TestStorage :
        
    """
    mongo DB scheme:
    
    
    operations:
    add_branch - look for given branch name, get and return last commit or create the new one
    append_to_branch - check again if we have new branch of same name ( first commit in the list do not match last one), create new branch or append
    (possible problem with preallocation namely for master branch)
    get_branches - select all branches, sorted by date of their last commit; find every commit 
    add_host - find or add config, find or add host; add possibly new config to host
    get_hosts - select all hosts; insert configs into them
    get_whole_tree(commit) - select all test cases of given commit_id; for every test case select all test runs and insert data
    get_summary_data(node_list, commit_list, arch_list)
             - for every node in node_list and every in commit_list: select all test cases under "node", group their run_id; get appropriate test_runs
               project to commute scaling, group ...
               possibly do for all commits on the list and group by commit_id
               
    ... design in such a way that we can use unwind and group aggregations to get summary values
    
    
    Provides interface to whole database of test data, namely:
    
    Adding, modify data:
    * append new commit chain to given branch (append_to_branch)
    * get last stored commit of a branch (last_commot_of_branch_
    * add host
    
    Reading:
    * get hosts
    * get branches
    * get whole test data  tree for given commit hash
    

Dashboard:
* get cumulative data for given nodes/blocks and given arch selection and every commit in given range
  cumulative data: for every given block and commit : totals for individual statuses (number of runs with such result 
  possibly filter benchmarks 
      data members:
      mc - mongo client
      
    """
    
    # Possible status values
    COMPILE_FAIL =20
    LINK_FAIL = 19
    RUN_FAIL = 18
    OUTPUT_FAIL = 17
    PERF_TIME_FAIL = 10
    PERF_MEM_FAIL = 9
    PERF_CALLS_FAIL =8
    PERF_SCALE_FAIL =7
    NO_FAIL = 0
    
    def __init__(self):
        self.mc = pymongo.MongoClient('localhost', 27017)
        self.db = self.mc['test_data']
        """
        branch = { name, remote, last_commit_hash, last_commit_date,  commits : [ commit_id, ... ], 
          # optional branch point
          root : { branch_id, commit_id }  }
        """
        self.branch = self.db['branch']    
        """
        commit = { id_, date, hash, author, message,
          # optional
          exam : { exam_id ,  branch_id } }      # in which branch we perform the exam
        """
        self.commit = self.db['commit']  
        """
        host = { os, hostname, uname, frequency, benchmark,
                config : [ config_id, .. ]}
        }
        """
        self.host = self.db['host']
        """
        config = { name, file_hash, file }
        """
        self.config = self.db['config']
        """
        test_case = { node_path, block_path,  commit_id
          runs: [ test_run_id, ..] }
        """
        self.test_case = self.db['test_case']
        """
        test_run = { host_id, config_id, commit_id,
          test_case_id, nproc, size, run_id,              
          time, memory, calls, status }
        """
        self.test_run = self.db['test_run']
        
      
    def add_branch(self,branch_info) :
        """
        Possibly add new branch (includes branch name, remote, and possibly address of web interface)
        For branch_info see gtest_parser.py, function get_branch_info.
        Returns last stored commit or None.
        """
        print "** add_branch\n", json.dumps(branch_info)
        result = self.branch.find( spec= { "name": branch_info['git_branch'] },
                                       fields = { "commits" : False } )
        if ( result.count() == 0 ) :
            last_hash = branch_info['git_commit']['hash']
            last_date = dateutil.parser.parse( branch_info['git_commit']['date_str'] )
            # make new branch
            self.branch.insert( { 'name' : branch_info['git_branch'],
                                  'remote' : branch_info['git_remote_branch'],
                                  'last_commit_hash' :last_hash,
                                  'last_commit_date' :last_date,
                                  'commits' : []
                                 } ) 
        elif ( result.count() ==1 ) :
            for branch in result:
                return branch['last_commit_hash']
        else :
            raise(ExcAssert("Non-unique branch name."))
      
    def append_to_branch(self,branch_info, commits):
        """
        Append array of commits to given branch. One commit is dict with structure:
        { 
          hash : git hash
          date : date and time of commit, commits are sorted by date 
          author : commit;s author
          message : whole commit message
          parent : { branch, commit_id } - very first commit of branch may have reference to origin of the branch
          This way commits are duplicated, but usually only two times one time in feature branch and second time in maser
          after merge.
        }  
        """
        print "** append_to_branch", json.dumps(branch_info), json.dumps(commits)
        pass
      
    def get_branches(self,depth=0):
        """
        Returns list of all branches with their commits. Optional parameter gives maximum 
        length of commit lists.
        """
        pass
      
    
    def add_host(self,host_data):
        """
        Add arch to existing host or create a new one.
        Returns error message if specification differs.
        host_data = 
        {
          os:               # {win|linux}_{32|64}
          hostname:         # local host name, address 
          uname:
          frequency:        # in MHz
          benchmark:        # time of sample code snipet
          arch : { name, file_hash, file }   # architecture - cmake.config 
        """
        print "** add_host\n",json.dumps(host_data)
        pass

    def get_hosts(self):
        """ 
        Returns list of all hosts with their architectures
        """
        pass
      
    def add_report(self,host_info,branch_info, exam_data) :
        self.add_host(host_info)
        print "** add_report\n", json.dumps(exam_data) 
        pass

    def get_whole_tree(self,commit):
        """
        Return dictionary with whole data tree of given commit.
        """
        pass
      
    def get_summary_data(self,node_list, commit_list, arch_list):
        """
        Return list of lists (first index in node_list, second in commit_list) with summary metrics over
        particular node childs and over arch_list. 
        one summary:
        { 
          status : { compile, link, run, output, perf, OK }
          time :
          memory :
          scale : { 8, 16 }
          calls :
        """
        pass
      


class BenchRequestHandler(BaseHTTPServer.BaseHTTPRequestHandler):
    
    
    
    def do_POST(self):
        srv=self.server
        response={}
        response_code=200
        content_type = self.headers['content-type']
        if ( content_type == 'application/json') :
            content_len = int(self.headers['content-length'])
            content = json.loads(self.rfile.read(content_len))
            if (content['request-type'] == 'add_report') :
                srv.context.host_info=content['host_info']
                srv.context.branch_info=content['branch_info']
                srv.context.exam_data = content['exam_data']
                last_commit=srv.ts.add_branch(srv.context.branch_info)
                response={ 'last_commit':last_commit }
            
            elif (content['request-type'] == 'add_to_branch') :
                srv.ts.append_to_branch(srv.context.branch_info,content['branch_data'])
                srv.ts.add_report(self.server.context.host_info, 
                                   self.server.context.branch_info,
                                   self.server.context.exam_data)                
            
            elif (content['request-type'] == 'get_branches') :
                response=srv.ts.get_branches( content['commit_list_length'] )
            
            elif (content['request-type'] == 'get_hosts') :    
                response=srv.ts.get_hosts()
                
            elif (content['request-type'] == 'get_node_tree') :    
                response=srv.ts.get_whole_tree(content['commit_id'])
               
            elif (content['request-type'] == 'get_summary') :    
                response=srv.ts.get_summary_data(content['node_id_list'],
                                                 content['commit_id_list'],
                                                 content['arch_id_list'])
            else :
                response_code  = 415
        else :
            response_code = 415

        self.send_response(response_code)
        self.send_header("Content-type", "application/json")
        self.end_headers()
        self.wfile.write(json.dumps(response))
        

class BenchHTTPServer(BaseHTTPServer.HTTPServer) :
    class context:
        host_info={}
        branch_info={}
        exam_data={}
        
    def __init__(self, address, storage) :
        BaseHTTPServer.HTTPServer.__init__(self,address, BenchRequestHandler)
        self.ts=storage
        
       
      
      
def start_server():
    server_address = ('', 8000)
    storage=TestStorage()
    httpd = BenchHTTPServer(server_address, storage)
    httpd.serve_forever() 
    
    
# Main
start_server()