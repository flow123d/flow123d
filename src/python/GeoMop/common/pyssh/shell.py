# -*- coding: utf-8 -*-

from __future__ import unicode_literals

import warnings
import ctypes

from . import api
from . import compat
from . import exceptions as exp


class Shell(object):
    """
    Shell session.
    """

    _channel = None

    def __init__(self, session, pty_size, env):
        self.session_wrapper = session
        self.session = session.session
        self.pty_size = pty_size
        self.env = env

    @property
    def channel(self):
        if self._channel is not None:
            return self._channel

        self._channel = api.library.ssh_channel_new(self.session);

        # Open ssh session
        ret = api.library.ssh_channel_open_session(self._channel)
        if ret != api.SSH_OK:
            raise exp.ConnectionError("Error code: {0}".format(ret))

        # Request pty
        ret = api.library.ssh_channel_request_pty(self._channel)
        if ret != api.SSH_OK:
            raise exp.ConnectionError("Error code: {0}".format(ret))

        # Request pty size
        #ret = api.library.ssh_channel_request_pty_size(self._channel,
        #                            self.pty_size[0], self.pty_size[1])
        #if ret != api.SSH_OK:
        #    raise RuntimeError("Error code: {0}".format(ret))

        # Request shell
        ret = api.library.ssh_channel_request_shell(self._channel)
        if ret != api.SSH_OK:
            raise exp.ConnectionError("Error code: {0}".format(ret))

        # Set environ variable if theese are available
        if self.env:
            for key, value in self.env.items():
                _key, _value = key, value
                if isinstance(_key, compat.text_type):
                    _key = compat.to_bytes(_key, encoding="utf-8")

                if isinstance(_value, compat.text_type):
                    _value = compat.to_bytes(_value, encoding="utf-8")

                res = api.library.ssh_channel_request_env(self.channel, _key, _value)
                res = api.library.ssh_channel_request_shell(self.channel)
                if res != api.SSH_OK:
                    msg = api.library.ssh_get_error(self.session)
                    print("Error: ", msg)
                    warnings.warn("Error on set {0} variable".format(key), RuntimeWarning)

        return self._channel

    def write(self, data):
        """
        Write bytes to remote shell.

        The `data` parameter accept both str and bytes, if you passes str (unicode) is
        automatically converted to bytes using utf-8 encoding.

        :param bytes data: arbitrary length of bytes.
        :returns: a number of bytes written to remote shell.
        :rtype: int
        """
        if isinstance(data, compat.text_type):
            data = compat.to_bytes(data, "utf-8")

        written = api.library.ssh_channel_write(self.channel, data, len(data))
        if written != len(data):
            raise RuntimeError("Error on write")
        return written

    def read(self, n):
        """
        Read bytes from remote shell.

        This method always return value, if not bytes available to read
        it returns an empty bytestring.

        :param int n: number of bytes to read
        :returns: bytestring of readed data.
        :rtype: bytes
        """

        res = api.library.ssh_channel_is_open(self.channel)
        if res == 0:
            raise RuntimeError("Channel is closed")

        res = api.library.ssh_channel_is_eof(self.channel)
        if res != 0:
            return b""

        buffer = ctypes.create_string_buffer(n)
        readed = api.library.ssh_channel_read_nonblocking(self.channel, buffer, n, 0)
        if readed < 0:
            raise RuntimeError("Error on read")

        return buffer.value

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        if self._channel is not None:
            if api.library.ssh_channel_is_closed(self._channel) == 0:
                api.library.ssh_channel_close(self._channel)
                api.library.ssh_channel_send_eof(self._channel)

            api.library.ssh_channel_free(self.channel)
            self._channel = None
