# -*- coding: utf-8 -*-

from __future__ import unicode_literals
import warnings

# def read_until_fail(sock):
#     retval = []
#     more_data=True
#     while more_data:
#         try:
#             read_data = sock.recv(1024)
#             retval.append(read_data)
#             more_data = len(read_data)>0
#         except:
#             more_data=False
#     return "".join(retval)


#class Forward(object):
#    channel = None
#
#    def __init__(self, session, local_host, local_port,  remote_host, remote_port):
#        self.local_port = local_port
#        self.local_host = local_host
#        self.remote_port = remote_port
#        self.remote_host = remote_host
#        self.session_wrapper = session
#        self.session = session.session
#
#    def start(self):
#        self.channel = api.library.ssh_channel_new(self.session)
#        if self.channel is None:
#            raise RuntimeError("Error open new channel")
#
#        ret = api.library.ssh_channel_open_forward(self.channel,
#            self.remote_host, self.remote_port, self.local_host, self.local_port)
#
#        if ret != api.SSH_OK:
#            raise RuntimeError("Error code: {0}".format(ret))
#
#
#        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
#        sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
#        sock.bind((self.local_host, self.local_port))
#        sock.listen(100)
#
#    def stop(self):
#        if self.channel is not None:
#            api.library.ssh_channel_close(self.channel)
#            api.library.ssh_channel_free(self.channel)
#            self.channel = None
#
#    def __del__(self):
#        self.stop()
