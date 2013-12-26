import BaseHTTPServer

class BenchRequestHandler(BaseHTTPServer.BaseHTTPRequestHandler):
    def do_POST(self):
        print "Header -----------"
        print self.headers
        content_len = int(self.headers.getheader('content-length'))
        data = self.rfile.read(content_len)
        
        print "Data -----------"
        print json.loads(data)

        self.send_response(201)
        self.end_headers()
        self.wfile.write("Report accepted.")
        
        
        
      
      
      
def start_server():
    server_address = ('', 8000)
    httpd = BaseHTTPServer.HTTPServer(server_address, BenchRequestHandler)
    httpd.serve_forever() 
    
    
# Main
start_server()