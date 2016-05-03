# -*- coding: utf-8 -*-

from __future__ import unicode_literals
import ctypes
import sys

from . import api
from . import compat


class LazyResult(object):
    """
    Lazy command execution result wrapper.

    This wrapper implements a iterator interface.
    """

    _return_code = None
    _consumed = False

    def __init__(self, session, command):
        self.session = session
        self.command = command

    def __next__(self):
        if self._finished:
            raise StopIteration()

        data = ctypes.create_string_buffer(10)
        readed_bytes = api.library.ssh_channel_read(self.channel, ctypes.byref(data),
                                                    len(data), 0)
        if readed_bytes > 0:
            return data.value

        api.library.ssh_channel_send_eof(self.channel);
        self._return_code = api.library.ssh_channel_get_exit_status(self.channel)
        api.library.ssh_channel_free(self.channel)
        self.channel = None
        self._finished = True
        raise StopIteration

    if sys.version_info[0] == 2:
        next = __next__

    def __iter__(self):
        if self._consumed:
            raise RuntimeError("Result are consumed")

        self._consumed = True
        self._finished = False

        self.channel = api.library.ssh_channel_new(self.session);

        # Open ssh session
        ret = api.library.ssh_channel_open_session(self.channel)
        if ret != api.SSH_OK:
            raise RuntimeError("Error code: {0}".format(ret))

        # Execute the command
        ret = api.library.ssh_channel_request_exec(self.channel, self.command)
        if ret != api.SSH_OK:
            msg = api.library.ssh_get_error(self.session)
            raise RuntimeError("Error {0}: {1}".format(ret, msg.decode('utf-8')))

        return self

    def as_bytes(self):
        """
        Launch the command and return a result as bytes.

        :returns: bytes chunk of command execution result
        :rtype: bytes
        """

        return b"".join([x for x in self])

    def as_str(self):
        """
        Launch the command and return a result as unicode string

        :returns: unicode chunk of command execution result
        :rtype: str/unicode
        """
        return self.as_bytes().decode("utf-8")

    def wait(self):
        """
        Waits a complete command execution and returns the return code

        :returns: execution result return code
        :rtype: int
        """
        list(self)
        return self.return_code

    @property
    def return_code(self):
        return self._return_code


class Result(LazyResult):
    """
    Consumed version of LazyResult. Useful for simple command
    execution.
    """
    _data = None

    def __init__(self, *args, **kwargs):
        super(Result, self).__init__(*args, **kwargs)

        # consume iterator and save state
        self._data = list(self)

    def as_bytes(self):
        """
        Return a cached result.

        :returns: bytes chunk of command execution result
        :rtype: bytes
        """
        return b"".join(self._data)

    def wait(self):
        return self.return_code
