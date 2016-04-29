# -*- coding: utf-8 -*-

from __future__ import unicode_literals
import warnings

from .session import Session
from .sftp import Sftp


def new_session(hostname="localhost", port="22", username=None,
                password=None, passphrase=None, connect_on_init=False,
                verify_knownhost_callback=None):
    """
    Shortcut method for create new session instance.

    Session by default has lazy connection management. It only connects
    when it is need. But this this function you can pass `connect_on_init`
    parameter with True and session connects to the remote server before
    it are returned.

    :param str hostname: remote ip or host
    :param int port: remote port
    :param str username: remote user name with which you want to authenticate
    :param str password: remote user password.
    :param str passphrase: passphrase in case you would authenticate with pubkey
    :param bool connect_on_init: determines the lazyness of connection with
                                 remote server.
    """

    session = Session(hostname=hostname, port=port, username=username,
                      password=password, passphrase=passphrase,
                      verify_knownhost_callback=verify_knownhost_callback)
    if connect_on_init:
        session._connect_if_not_connected()
    return session


def connect(hostname="localhost", port="22", username=None,
            password=None, passphrase=None):
    """
    Shortcut method for create new session and connects to remote server.

    :param str hostname: remote ip or host
    :param int port: remote port
    :param str username: remote user name with which you want to authenticate
    :param str password: remote user password.
    :param str passphrase: passphrase in case you would authenticate with pubkey

    **NOTE:** this method is deprecated.
    """

    warnings.warn("connect function is deprectad, use new_session function",
                  DeprecationWarning)
    return new_session(hostname=hostname, port=port, username=username,
                       password=password, passphrase=passphrase, connect_on_init=True)
