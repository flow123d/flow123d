#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
# ----------------------------------------------
import importlib
import json
import platform
# ----------------------------------------------
import yaml

from scripts.core.base import Paths, Printer
# ----------------------------------------------


job_ok_string = '[JOB_STATUS:OK]'


class DummyModule(object):
    """
    Dummy class for auto-completion purposes in IDE
    :type template: str
    """

    def __init__(self):
        self.template = None

    def Module(self, test_case, proc_value, filename):
        """
        """
        pass

    def ModuleJob(self, job_id):
        """
        :rtype : scripts.pbs.job.Job
        """
        pass


def get_pbs_module(hostname_hint=None):
    """
    file host_table.yaml serves as lookup table when using python script in queue mode
    each key is hostname and each value names a module which should be loaded
    modules are located in /src/python/scripts/pbs/modules

    If no matching key for current machine exists try to use pbs_<hostname>
    where all dots(.) are replaced with underscores(_)

    if hostname_hint is not set node name will be used
    :rtype : scripts.pbs.modules.pbs_tarkil_cesnet_cz
    """
    pbs_module_path = None
    host_file = Paths.join(Paths.flow123d_root(), 'config', 'host_table.yaml')
    host_file_exists = Paths.exists(host_file)
    hostname = hostname_hint or platform.node()
    from_host = False

    # try to get name from json file
    if host_file_exists:
        with open(host_file, 'r') as fp:
            hosts = yaml.load(fp)
            pbs_module_path = hosts.get(hostname, None)
            from_host = pbs_module_path is not None

    if not pbs_module_path:
        hostname = hostname.replace('.', '_')
        pbs_module_path = 'pbs_{}'.format(hostname)

    # construct full path for import
    full_module_path = 'scripts.pbs.modules.{module_name}'.format(module_name=pbs_module_path)

    # try to get pbs_module
    try:
        return importlib.import_module(full_module_path)
    except ImportError:
        Printer.err('Error:  Could not load module "{}" ({}) for hostname "{}"',
                    pbs_module_path, full_module_path, hostname)
        Printer.open(2)
        if host_file_exists:
            if from_host:
                Printer.err('Value specified in host_table.yaml "{}" points to non-existing module', pbs_module_path)
            else:
                Printer.err('Config file host_table.yaml does not have entry for hostname "{}"', hostname)
        else:
            Printer.err('Config file host_table.yaml does not exists ({}) and auto module detection failed', host_file)
        Printer.close(2)
        raise
