#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

"""
This script is simple iface to coverage_merge_module module.
Script is expected to run as main python file. Script takes arguments from console and calls methods of module.

Usage: coverage_merge_script.py [options] [file1 file2 ... filen]

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -o FILE, --output=FILE
                        output file xml name
  -p FILE, --path=FILE  xml location, default current directory
  -l LOGLEVEL, --log=LOGLEVEL
                        Log level DEBUG, INFO, WARNING, ERROR, CRITICAL
  -f, --filteronly      If set all files will be filtered by keep rules
                        otherwise all given files will be merged and filtered.
  -s SUFFIX, --suffix=SUFFIX
                        Additional suffix which will be added to filtered
                        files so they original files can be preserved
  -k NAME, --keep=NAME  preserves only specific packages. e.g.:
                        'python merge.py -k src.la.*'
                        will keep all packgages in folder src/la/ and all
                        subfolders of this folders.
                        There can be mutiple rules e.g.:
                        'python merge.py -k src.la.* -k unit_tests.la.'
                        Format of the rule is simple dot (.) separated names
                        with wildcard (*) allowed, e.g:
                        package.subpackage.*

If no files are specified all xml files in current directory will be selected.
Useful when there is not known precise file name only location
"""

from __future__ import absolute_import

import pathfix
pathfix.append_to_path()

import system.versions
system.versions.require_version_2()

from optparse import OptionParser
from coverage import coverage_merge_module


# parse arguments
def create_parser ():
    """Creates command line parse"""
    newline = 10 * '\t';
    parser = OptionParser (usage="%prog [options] [file1 file2 ... filen]", version="%prog 1.0",
                           epilog="If no files are specified all xml files in current directory will be selected. \n" +
                                  "Useful when there is not known precise file name only location")

    parser.add_option ("-o", "--output", dest="filename", default="coverage-merged.xml",
                       help="output file xml name", metavar="FILE")
    parser.add_option ("-p", "--path", dest="path", default="./",
                       help="xml location, default current directory", metavar="FILE")
    parser.add_option ("-l", "--log", dest="loglevel", default="DEBUG",
                       help="Log level DEBUG, INFO, WARNING, ERROR, CRITICAL")
    parser.add_option ("-f", "--filteronly", dest="filteronly", default=False, action='store_true',
                       help="If set all files will be filtered by keep rules otherwise " +
                            "all given files will be merged and filtered.")
    parser.add_option ("-s", "--suffix", dest="suffix", default='',
                       help="Additional suffix which will be added to filtered files so they original files can be preserved")
    parser.add_option ("-k", "--keep", dest="packagefilters", default=None, metavar="NAME", action="append",
                       help="preserves only specific packages. e.g.: " + newline +
                            "'python merge.py -k src.la.*'" + newline +
                            "will keep all packgages in folder " +
                            "src/la/ and all subfolders of this folders. " + newline +
                            "There can be mutiple rules e.g.:" + newline +
                            "'python merge.py -k src.la.* -k unit_tests.la.'" + newline +
                            "Format of the rule is simple dot (.) separated names with wildcard (*) allowed, e.g: " + newline +
                            "package.subpackage.*")
    return parser


def parse_args (parser):
    """Parses argument using given parses and check resulting value combination"""
    (options, args) = parser.parse_args ()

    # for now, no check needed

    return (options, args)


if __name__ == '__main__':
    parser = create_parser ()
    (options, args) = parse_args (parser)

    coverage_merge_module.merge (options, args)

