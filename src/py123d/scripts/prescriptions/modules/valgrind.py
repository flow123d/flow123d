#!/usr/bin/python
# -*- coding: utf-8 -*-
# author: Jan Hybs

from py123d.loggers import printf
from py123d.scripts.core.execution import BinExecutor, OutputMode


def massif_hook(pypy):
    """
    Function will initialize ValgrindMassif and print result on success
    :type pypy: scripts.core.pypy.PyPy
    """
    vm = ValgrindMassif(pypy.case.fs.valgrind_out)
    vm.output = OutputMode.variable_output()
    vm.start()
    vm.join()

    # if ms_print ended successfully we print first 32 lines (where graph is)
    if vm.returncode == 0:
        output = vm.output.read()
        lines = str(output).splitlines()[0:32]
        printf.out('\n    '.join(lines))


class ValgrindMassif(BinExecutor):
    """
    Class ValgrindMassif is simple bin executor which will call ms_print and collect result
    """

    def __init__(self, massif_file):
        """
        :type massif_file: str
        """
        command = ['ms_print', massif_file]
        super(ValgrindMassif, self).__init__(command, 'massif-ms_print')
