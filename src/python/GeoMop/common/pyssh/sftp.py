# -*- coding: utf-8 -*-

from __future__ import unicode_literals

import ctypes
import stat
import os
import io

from . import api
from . import compat
from . import exceptions as exp


class Sftp(object):
    """
    Sftp wrapper.

    Exposes api for interacting with sftp subsystem: put or get files,
    open files with random read-write access, etc.

    :ivar ponter sftp: c sftp session pointer
    :ivar pointer session: c ssh session pointer

    :param pyssh.session.Session session: initialized and connected
        :py:class:`pyssh.session.Session` instance.
    """

    sftp = None
    session = None

    def __init__(self, session, buffer_size=1024*16):
        self.session_wrapper = session
        self.session = session.session

        self.buffer_size = buffer_size

        # TODO: handle exceptions
        self.sftp = api.library.sftp_new(self.session)
        api.library.sftp_init(self.sftp)

    def _get_file_metadata(self, file_ptr):
        attrs_ptr = api.library.sftp_fstat(file_ptr)
        if attrs_ptr is None:
            msg = api.library.ssh_get_error(self.session)
            raise exp.ConnectionError("Error raised by ssh: {0}".format(msg.decode("utf-8")))

        attrs = ctypes.cast(attrs_ptr, ctypes.POINTER(api.SftpAttributes))
        return attrs.contents

    def _open_remote_file(self, path):
        remote_file_ptr = api.library.sftp_open(self.sftp, path, os.O_RDONLY, stat.S_IRWXU)

        if remote_file_ptr is None:
            msg = api.library.ssh_get_error(self.session)
            raise exp.ConnectionError("Error raised by ssh: {0}".format(msg.decode("utf-8")))

        return remote_file_ptr

    def get(self, remote_path, local_path):
        """
        Get a remote file to local.

        :param str remote_path: remote file path
        :param str local_path:  local file path
        """
        remote_path = compat.to_bytes(remote_path)

        # Create new pointer to remote file
        remote_file_ptr = self._open_remote_file(remote_path)

        # Obtain metadata with remote file size
        remote_file_attrs = self._get_file_metadata(remote_file_ptr)

        stats = {
            "total_readed": 0,
            "total_size": remote_file_attrs.size,
            "errors_counter": 0,
        }

        def read_pipeline(f):
            while True:
                buffer = ctypes.create_string_buffer(self.buffer_size)
                readed = api.library.sftp_read(remote_file_ptr, ctypes.byref(buffer),
                                               self.buffer_size)
                if readed == 0:
                    if stats["total_readed"] != stats["total_size"]:
                        stats["errors_counter"] += 1
                        return False
                    return True

                elif readed < 0:
                    raise exp.ConnectionError("Connection interrumped")

                else:
                    stats["errors_counter"] = 0
                    stats["total_readed"] += readed
                    f.write(buffer.raw[:readed])

        try:
            with io.open(local_path, "wb") as f:
                while stats["errors_counter"] < 3:
                    ok = read_pipeline(f)
                    if ok:
                        break

                    api.library.sftp_close(remote_file_ptr)
                    remote_file_ptr = self._open_remote_file(remote_path)

                if stats["errors_counter"] >= 3:
                    raise exp.Connection("Connection errors repeated more than 3 times")

        except (exp.ConnectionError, RuntimeError):
            api.library.sftp_close(remote_file_ptr)
            raise

    def put(self, path, remote_path):
        """
        Puts the local file to remote host.

        :param str path: local file path
        :param str remote_path: remote file path
        """

        if not os.path.exists(path):
            raise RuntimeError("Path {0} does not exists".format(path))

        if isinstance(remote_path, compat.text_type):
            remote_path = compat.to_bytes(remote_path, "utf-8")

        access_type = os.O_WRONLY | os.O_CREAT | os.O_TRUNC
        remote_file_ptr = api.library.sftp_open(self.sftp, remote_path, access_type, stat.S_IRWXU)

        if remote_file_ptr is None:
            msg = api.library.ssh_get_error(self.session)
            raise exp.ConnectionError("Error raised by ssh: {0}".format(msg.decode("utf-8")))

        with io.open(path, "rb") as f:
            while True:
                chuck = f.read(self.buffer_size)
                if not chuck:
                    break

                written = api.library.sftp_write(remote_file_ptr, chuck, len(chuck))
                if written != len(chuck):
                    raise RuntimeError("Can't write file")

        api.library.sftp_close(remote_file_ptr)

    def open(self, path, mode):
        """
        Open a remote file.

        :param str path: remote file path
        :param int mode: open file model
                         (see http://docs.python.org/3.3/library/os.html#open-flag-constants)

        :returns: SFTP File wrapper
        :rtype: pyssh.SftpFile
        """
        if isinstance(path, compat.text_type):
            path = compat.to_bytes(path, "utf-8")

        return SftpFile(path, mode, self)

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        api.library.sftp_free(self.sftp)


class SftpFile(object):
    """
    SFTP File wrapper
    """

    _closed = False

    def __init__(self, path, mode, sftp_wrapper):
        self.sftp_wrapper = sftp_wrapper
        self.sftp = sftp_wrapper.sftp

        self.file = api.library.sftp_open(self.sftp, path, mode, stat.S_IRWXU)

        if self.file is None:
            self._closed = True
            raise ext.SftpError("Can't open file {0}".format(path.decode("utf-8")))

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        self.close()

    def write(self, data):
        """
        Write bytes to remote file.

        :param bytes data: bytes chunk of data
        :returns: number of bytes are written
        :rtype: int
        """
        written = api.library.sftp_write(self.file, data, len(data))
        if written != len(data):
            raise RuntimeError("Can't write file")

        return written

    def read(self, num=None, buffer_length=1024):
        """
        Read from remote file.

        :param int num: number of bytes to read, if num is None reads all.
        :returns: readed bytes chunk
        :rtype: bytes
        """
        if num is None:
            buffer_len = buffer_length
        else:
            buffer_len = num

        buffer = ctypes.create_string_buffer(buffer_len)
        readed = api.library.sftp_read(self.file, ctypes.byref(buffer),  buffer_len);

        if readed == 0:
            return b""

        if num is not None and num > 0:
            if buffer_len != readed:
                raise RuntimeError("Error on read")
            return buffer.raw

        readed_data = [buffer.raw]
        while True:
            buffer = ctypes.create_string_buffer(buffer_len)

            readed = api.library.sftp_read(self.file, ctypes.byref(buffer),  buffer_len);
            if readed == 0:
                break

            readed_data.append(buffer.raw[:readed])
        return b"".join(readed_data)

    def seek(self, offset):
        """
        Change position on a remote file.

        :param int offset: file position
        :returns: boolean value if seek is success or not
        :rtype: bool
        """
        ret = api.library.sftp_seek64(self.file, offset);
        if ret != api.SSH_OK:
            return False

        return True

    def tell(self):
        """
        Query the current position on a file.

        :returns: a current position.
        :rtype: int
        """
        return api.library.sftp_tell64(self.file)

    def close(self):
        """
        Close a opened file.
        """
        if self._closed:
            raise ext.ResourceManagementError("SftpFile instance already closed.")

        self._closed = True
        api.library.sftp_close(self.file)
