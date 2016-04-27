# -*- coding: utf-8 -*-

from __future__ import unicode_literals
import functools
import warnings

from . import api
from . import compat
from . import result
from . import exceptions as exp

from . import shell
from . import sftp
from ctypes import byref, c_char_p


def _lazy_connect(func):
    @functools.wraps(func)
    def _decorator(self, *args, **kwargs):
        self._connect_if_not_connected()
        return func(self, *args, **kwargs)
    return _decorator


def _check_open_session(func):
    @functools.wraps(func)
    def _decorator(self, *args, **kwargs):
        if self._closed:
            raise exp.SshError("Session aleady closed.")
        return func(self, *args, **kwargs)
    return _decorator


class Session(object):
    """
    SSH Session wrapper.

    Actually accepts two methods for authentication: the simple a simple password or
    a pubkey. If password is not provided, attempts using pubkey, with or without pasphrase.

    :ivar pointer session: c ssh session pointer
    :ivar bytes username: current username

    :param str hostname: remote ip or host
    :param int port: remote port
    :param str username: remote user name with which you want to authenticate
    :param str password: remote user password.
    :param str passphrase: passphrase in case you would authenticate with pubkey
    :param func verify_knownhost_callback: function which gets called upon connecting to host. Should return
        True if connection is allowed, False otherwise. The only parameter to the function is remote host key 
        SHA1 hash. WARNING: you should always verify host signature!
    """

    session = None
    username = None
    password = None

    _closed = False
    _connected = False

    def __init__(self, hostname, port=22, username=None, password=None, passphrase=None, verify_knownhost_callback=None):
        self.session = api.library.ssh_new()

        if isinstance(hostname, compat.text_type):
            self.hostname = compat.to_bytes(hostname)
        else:
            self.hostname = hostname

        if isinstance(port, int):
            self.port = compat.to_bytes(str(port))
        elif isinstance(port, compat.text_type):
            self.port = compat.to_bytes(port)
        else:
            self.port = port

        if isinstance(username, compat.text_type):
            self.username = compat.to_bytes(username)
        else:
            self.username = username

        if isinstance(password, compat.text_type):
            self.password = compat.to_bytes(password)
        else:
            self.password = password

        if self.username:
            api.library.ssh_options_set(self.session, api.SSH_OPTIONS_USER, self.username)

        if isinstance(passphrase, compat.text_type):
            self.passphrase = compat.to_bytes(passphrase, "utf-8")
        else:
            self.passphrase = passphrase

        api.library.ssh_options_set(self.session, api.SSH_OPTIONS_PORT_STR, self.port)
        api.library.ssh_options_set(self.session, api.SSH_OPTIONS_HOST, self.hostname)

        self.verify_knownhost_callback = verify_knownhost_callback

    def _connect_if_not_connected(self):
        # Do nothing if it is connected
        if self._connected:
            return

        self._connected = True

        ret = api.library.ssh_connect(self.session)
        if ret != api.SSH_OK:
            remote_msg = compat.to_text(api.library.ssh_get_error(self.session))
            msg = ("Unable to connect to remote server. "
                   "(Return code: {0}, Return message: {1})")

            raise exp.ConnectionError(msg.format(ret, remote_msg))

        if self.verify_knownhost_callback is not None:
            hash = c_char_p()
            try:
                hashlen = api.library.ssh_get_pubkey_hash(self.session, byref(hash))
                if hashlen < 1:
                    raise exp.HostVerificationError("Error verifying remote host - could not fetch pubkey hash.")
                if not self.verify_knownhost_callback(hash.value[0:hashlen]):
                    raise exp.HostVerificationError("Error verifying remote host - host not authentic.")
            finally:
                api.library.ssh_clean_pubkey_hash(hash)

        if self.password is not None:
            ret = api.library.ssh_userauth_password(self.session, None, self.password)
            if ret != api.SSH_AUTH_SUCCESS:
                raise exp.AuthenticationError("Error when trying authenticate with password. "
                                              "(Error code: {0})".format(ret))
        else:
            ret = api.library.ssh_userauth_autopubkey(self.session, self.passphrase)
            if ret != api.SSH_AUTH_SUCCESS:
                raise exp.AuthenticationError("Error when trying authenticate with pubkey. "
                                              "(Error code: {0})".format(ret))
    def close(self):
        """
        Close initialized ssh connection.
        """
        if self._closed:
            raise exp.ResourceManagementError("Already closed")

        self._closed = True

        if self._connected:
            api.library.ssh_disconnect(self.session)
            self._connected = False

        api.library.ssh_free(self.session)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    @_check_open_session
    @_lazy_connect
    def create_shell(self, pty_size=(80, 24), env={}):
        """
        Creates a new shell session  throught current ssh channel.

        :param tuple pty_size: in case of shell is true this indicates
            the size of a virtual terminal
        :param dict env: addiotional environ variables
        """
        # warnings.warn("Shell feature is very experimental and uncomplete.", Warning)
        return shell.Shell(self, pty_size, env)

    @_check_open_session
    @_lazy_connect
    def create_sftp(self):
        """
        Create a new sftp session throught current ssh channel.

        :returns: Sftp instance
        :rtype: :py:class:`pyssh.sftp.Sftp`
        """
        return sftp.Sftp(self)

    @_check_open_session
    @_lazy_connect
    def execute(self, command, lazy=False):
        """
        Execute command on remote host.

        This command can return :py:class:`~pyssh.result.Result` or
        :py:class:`~pyssh.result.LazyResult` depending of lazy parameter.

        :param str command: command string
        :param bool lazy: set true for return a lazy result
                          instead a evaluated. Useful for execute
                          commands with large output (default: False)

        :returns: Result instance
        :rtype: :py:class:`pyssh.result.Result`
        """

        if isinstance(command, compat.text_type):
            command = compat.to_bytes(command)

        if lazy:
            _result = result.LazyResult(self.session, command)
        else:
            _result = result.Result(self.session, command)
        return _result
