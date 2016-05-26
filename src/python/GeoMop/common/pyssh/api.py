# -*- coding: utf-8 -*-

import warnings
import ctypes
import ctypes.util


def load_library():
    libpath = ctypes.util.find_library('ssh')
    libssh = ctypes.CDLL(libpath)
    return libssh


SSH_OK = 0
SSH_ERROR = -1
SSH_AGAIN = -2
SSH_EOF = -127

SSH_OPTIONS_HOST = 0
SSH_OPTIONS_PORT = 1
SSH_OPTIONS_PORT_STR = 2
SSH_OPTIONS_FD = 3
SSH_OPTIONS_USER = 4
SSH_OPTIONS_SSH_DIR = 5
SSH_OPTIONS_IDENTITY = 6
# TODO...

SSH_AUTH_SUCCESS = 0
SSH_AUTH_DENIED = 1
SSH_AUTH_PARTIAL = 2
SSH_AUTH_INFO = 3
SSH_AUTH_AGAIN = 4
SSH_AUTH_ERROR = -1


class SftpAttributes(ctypes.Structure):
    _fields_ = [("name", ctypes.c_char_p),
                ("longname", ctypes.c_char_p),
                ("flags", ctypes.c_uint32),
                ("type", ctypes.c_uint8),
                ("size", ctypes.c_uint64),]


try:
    library = load_library()
    library.ssh_new.argtypes = []
    library.ssh_new.restype = ctypes.c_void_p

    library.ssh_free.argtypes = [ctypes.c_void_p]

    library.ssh_connect.argtypes = [ctypes.c_void_p]
    library.ssh_connect.restype = ctypes.c_int

    library.ssh_disconnect.argtypes = [ctypes.c_void_p]
    library.ssh_options_set.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p]

    library.ssh_userauth_password.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_char_p]
    library.ssh_userauth_password.restype = ctypes.c_int

    library.ssh_userauth_autopubkey.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
    library.ssh_userauth_autopubkey.restype = ctypes.c_int

    library.ssh_channel_new.argtypes = [ctypes.c_void_p]
    library.ssh_channel_new.restype = ctypes.c_void_p

    library.ssh_channel_open_session.argtypes = [ctypes.c_void_p]
    library.ssh_channel_open_session.restype = ctypes.c_int

    library.ssh_channel_request_exec.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
    library.ssh_channel_request_exec.restype = ctypes.c_int

    library.ssh_channel_read.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_uint, ctypes.c_int]
    library.ssh_channel_read.restype = ctypes.c_int

    library.ssh_channel_read_nonblocking.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_uint, ctypes.c_int]
    library.ssh_channel_read_nonblocking.restype = ctypes.c_int

    library.ssh_channel_write.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_uint]
    library.ssh_channel_write.restype = ctypes.c_int

    library.ssh_channel_send_eof.argtypes = [ctypes.c_void_p]
    library.ssh_channel_send_eof.restype = ctypes.c_int

    library.ssh_channel_is_eof.argtypes = [ctypes.c_void_p]
    library.ssh_channel_is_eof.restype = ctypes.c_int

    library.ssh_channel_is_open.argtypes = [ctypes.c_void_p]
    library.ssh_channel_is_open.restype = ctypes.c_int

    library.ssh_channel_is_closed.argtypes = [ctypes.c_void_p]
    library.ssh_channel_is_closed.restype = ctypes.c_int

    library.ssh_channel_close.argtypes = [ctypes.c_void_p]
    library.ssh_channel_close.restype = ctypes.c_int

    library.ssh_channel_free.argtypes = [ctypes.c_void_p]

    library.ssh_channel_get_exit_status.argtypes = [ctypes.c_void_p]
    library.ssh_channel_get_exit_status.restype = ctypes.c_int

    library.ssh_channel_request_env.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_char_p]
    library.ssh_channel_request_env.restype = ctypes.c_int

    library.ssh_channel_request_pty.argtypes = [ctypes.c_void_p]
    library.ssh_channel_request_pty.restype = ctypes.c_int

    library.ssh_channel_request_pty_size.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
    library.ssh_channel_request_pty_size.restype = ctypes.c_int

    library.ssh_channel_request_shell.argtypes = [ctypes.c_void_p]
    library.ssh_channel_request_shell.restype = ctypes.c_int

    library.ssh_get_error.argtypes = [ctypes.c_void_p]
    library.ssh_get_error.restype = ctypes.c_char_p

    library.ssh_get_pubkey_hash.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_char_p)]
    library.ssh_get_pubkey_hash.restype = ctypes.c_int

    library.ssh_clean_pubkey_hash.argtypes = [ctypes.POINTER(ctypes.c_char_p)]
    library.ssh_clean_pubkey_hash.restype = None

    # SFTP
    library.sftp_new.argtypes = [ctypes.c_void_p]
    library.sftp_new.restype = ctypes.c_void_p

    library.sftp_init.argtypes = [ctypes.c_void_p]
    library.sftp_init.restype = None

    library.sftp_free.argtypes = [ctypes.c_void_p]

    library.sftp_fstat.argtypes = [ctypes.c_void_p]
    library.sftp_fstat.restype = SftpAttributes
    library.sftp_fstat.restype = ctypes.c_void_p

    library.sftp_open.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int, ctypes.c_int]
    library.sftp_open.restype = ctypes.c_void_p

    library.sftp_close.argtypes = [ctypes.c_void_p]

    library.sftp_write.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_uint]
    library.sftp_write.restype = ctypes.c_int

    library.sftp_seek64.argtypes = [ctypes.c_void_p, ctypes.c_ulonglong]
    library.sftp_seek64.restype = ctypes.c_int

    library.sftp_tell64.argtypes = [ctypes.c_void_p]
    library.sftp_tell64.restype = ctypes.c_ulonglong

    library.sftp_read.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_uint]
    library.sftp_read.restype = ctypes.c_int

    # Forward
    library.ssh_channel_open_forward.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_int]
    library.ssh_channel_open_forward.restype = ctypes.c_int

except (AttributeError, OSError, IOError):
    warnings.warn("ssh shared library not found or incompatible")
