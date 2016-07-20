# -*- coding: utf-8 -*-

class SshError(Exception):
    pass


class SftpError(SshError):
    pass


class AuthenticationError(SshError):
    pass


class ConnectionError(SshError):
    pass


class ResourceManagementError(SshError):
    pass


class HostVerificationError(SshError):
    pass
