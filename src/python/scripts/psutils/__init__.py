#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

import platform
system = platform.system().lower()


if system.startswith('linux'):
    from scripts.psutils.linux_psutil import Process

# if system.startswith('windows'):
#     from scripts.psutils.windows_psutil import Process
#
# if system.startswith('cygwin'):
#     from scripts.psutils.cygwin_psutil import Process
