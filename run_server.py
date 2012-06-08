#/usr/bin/python

from CGIHTTPServer import CGIHTTPRequestHandler, BaseHTTPServer

if __name__ == '__main__' :

    HandlerClass = CGIHTTPRequestHandler
    ServerClass = BaseHTTPServer.HTTPServer
    HandlerClass.cgi_directories = ['/']

    #server_address = ('18.55.1.170',8000)
    server_address = ('0.0.0.0',8000)
    httpd = ServerClass(server_address,HandlerClass)

    sa = httpd.socket.getsockname()
    print "Serving HTTP on", sa[0], "port", sa[1], "..."
    httpd.serve_forever()
