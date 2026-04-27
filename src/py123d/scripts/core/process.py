#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
import time

import psutil


class ProcessUtils(object):
    KB = 1000.0
    MB = KB ** 2
    GB = KB ** 3

    KiB = 1024.0
    MiB = KiB ** 2
    GiB = KiB ** 3

    _reasonable_amount_of_time = 1
    _just_a_sec = 0.1

    @classmethod
    def list_children(cls, process):
        """
        :type process: psutil.Process
        """
        result = []
        for child in process.children():
            result.extend(cls.list_children(child))
        result.append(process)
        return result

    @classmethod
    def get_memory_info(cls, process, prop='vms', units=MiB):
        """
        :type process: psutil.Process
        """
        children = cls.list_children(process)
        total = 0
        for child in children:
            try:
                total += getattr(child.memory_info(), prop)
            except psutil.NoSuchProcess:
                continue
        return total / units

    @classmethod
    def apply(cls, children, prop, *args, **kwargs):
        result = []
        for child in children:
            try:
                result.append(getattr(child, prop)(*args, **kwargs))
            except psutil.NoSuchProcess as e: pass
        return result

    @classmethod
    def terminate(cls, process):
        cls.terminate_all(cls.list_children(process))

    @classmethod
    def terminate_all(cls, children):
        cls.apply(children, 'terminate')

    @classmethod
    def kill(cls, process):
        cls.kill_all(cls.list_children(process))

    @classmethod
    def kill_all(cls, children):
        cls.apply(children, 'kill')

    @classmethod
    def secure_kill(cls, process):
        # first, lets be reasonable and terminate all processes (SIGTERM)
        children = cls.list_children(process)
        cls.terminate_all(children)

        # wait jus a sec to let it all blend in
        time.sleep(cls._just_a_sec)

        # if some process are still running
        if True in cls.apply(children, 'is_running'):
            # wait a bit to they can safely exit...
            time.sleep(cls._reasonable_amount_of_time)
            # check status again, maybe perhaps they finish what
            # they have started and are done?
            if True in cls.apply(children, 'is_running'):
                # looks like they are still running to SIGKILL 'em
                cls.kill_all(children)
            else:
                # all processes finish they job in reasonable period since
                # SIGTERM was sent
                return True
        else:
            # all processes finish right after* SIGTERM was announced
            return True

        # return max value either True or False
        # False meaning some process is still running

        return not max(cls.apply(children, 'is_running'))

