#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

import enum
import os, datetime


class ICollectTool(object):
    """
    Class ICollectTool is abstract class for any collection tool
    """

    include = None
    exclude = None

    def process_file(self, f: str) -> list:
        raise NotImplemented('Method must be implemented')


class LogPolicy(enum.Enum):
    ON_ERROR = 1
    ALWAYS = 2
    NEVER = 3


class CollectResult(object):
    """
    Class CollectResult is simple crate for holding result from
     collect tool
    :type items           : list[dict]
    :type logs            : list[dict]
    :type logs_updated    : list[dict]
    """

    def __init__(self, items: list, logs: list, log_policy: LogPolicy=LogPolicy.ON_ERROR, log_folder: str=None):
        self.items = items
        self.logs = []
        self.logs_updated = []

        if log_policy is LogPolicy.NEVER:
            logs = []
        elif log_policy is LogPolicy.ALWAYS:
            pass
        elif log_policy is LogPolicy.ON_ERROR:
            try:
                rc = items[0].get('base', {}).get('returncode', 1)
            except:
                rc = 1

            if rc in (0, None):
                logs = []

        for log_file in logs:
            log = self._read_log(log_file)
            self.logs.append(dict(
                filename=log_file.replace(log_folder, '').strip('/'),
                data=log
            ))

    def update(self, log_ids: list, dest='logs'):
        """
        Method updates items with specific mongo ids
        :param log_ids:
        :param dest:
        :return:
        """
        self.logs_updated = []
        for i in range(len(log_ids)):
            self.logs_updated.append(dict(
                filename=self.logs[i].get('filename'),
                data=log_ids[i],
                filesize=len(self.logs[i].get('data')) if self.logs[i].get('data') else 0,
            ))

        for item in self.items:
            item[dest] = self.logs_updated

    @classmethod
    def _read_log(cls, log_file: str):
        if os.path.exists(log_file):
            return open(log_file, 'rb').read()
        return None
