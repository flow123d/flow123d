# encoding: utf-8
# author:   Jan Hybs

"""
This script is simple iface to ist_formatter_module module.
Script is expected to run as main python file. Script takes arguments from console and calls methods of module.

Usage: Usage: ist_script.py [options]

Options:
  -h, --help            show this help message and exit
  -i FILENAME, --input=FILENAME
                        Absolute or relative path to JSON file which will be
                        processed
  -o FILENAME, --output=FILENAME
                        Absolute or relative path output file which will be
                        generated/overwritten
  -f FORMAT, --format=FORMAT
                        'tex' or 'html' output

"""

import pathfix

pathfix.append_to_path()

import system.versions

system.versions.require_version_2()

import sys
from optparse import OptionParser
from ist.ist_formatter_module import ISTFormatter


def create_parser():
    """Creates command line parse"""
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="input", metavar="FILENAME", default=None,
                      help="Absolute or relative path to JSON file which will be processed")
    parser.add_option("-o", "--output", dest="output", metavar="FILENAME", default=None,
                      help="Absolute or relative path output file which will be generated/overwritten")
    parser.add_option("-f", "--format", dest="format", metavar="FORMAT", default="tex",
                      help="tex|html output format")
    parser.set_usage("""%prog [options]""")
    return parser


def parse_args(parser):
    """Parses argument using given parses and check resulting value combination"""
    (options, args) = parser.parse_args()

    if options.input is None:
        print "Error: No input file specified!"
        parser.print_help()
        sys.exit(1)

    if options.output is None:
        print "Error: No output file specified!"
        parser.print_help()
        sys.exit(1)

    return options, args


def main():
    parser = create_parser()
    (options, args) = parse_args(parser)

    # create instance of formatter
    formatter = ISTFormatter()

    # convert to tex format
    if options.format.lower() in ('tex', 'latex'):
        formatter.json2latex(options.input, options.output)
        sys.exit(0)

    # convert to HTML format
    if options.format.lower() in ('html', 'html5', 'www', 'htm'):
        formatter.json2html(options.input, options.output)
        sys.exit(0)

    print "Error: Unsupported format '{:s}'".format(options.format)
    sys.exit(1)


# only if this file is main python file, read args
if __name__ == "__main__":
    main()