# encoding: utf-8
# author:   Jan Hybs

"""
This script is simple iface to profiler_formatter_module module.
Script is expected to run as main python file. Script takes arguments from console and calls methods of module.

Usage: profiler_formatter_script.py [options]

Options:
  -h, --help            show this help message and exit
  -i FILENAME, --input=FILENAME
                        Absolute or relative path to JSON file which will be
                        processed
  -o FILENAME, --output=FILENAME
                        Absolute or relative path output file which will be
                        generated/overwritten
  -f CLASSNAME, --formatter=CLASSNAME
                        Classname of formatter which will be used, to list
                        available formatters use option -l (--list)
  -l, --list            Prints all formatters available in folder formatters
                        (using duck-typing)
  -s STYLES, --style=STYLES
                        Additional styling options in name:value format (for
                        example separator:  default is os separator)

"""

import pathfix

pathfix.append_to_path ()

import system.versions

system.versions.require_version_2 ()

import sys
from optparse import OptionParser
from profiler.profiler_formatter_module import ProfilerFormatter


def create_parser ():
    """Creates command line parse"""
    parser = OptionParser ()
    parser.add_option ("-i", "--input", dest="input", metavar="FILENAME", default=None,
                       help="Absolute or relative path to JSON file which will be processed")
    parser.add_option ("-o", "--output", dest="output", metavar="FILENAME", default=None,
                       help="Absolute or relative path output file which will be generated/overwritten")
    parser.add_option ("-f", "--formatter", dest="formatter", metavar="CLASSNAME", default="SimpleTableFormatter",
                       help="Classname of formatter which will be used, to list available formatters use option -l (--list)")
    parser.add_option ("-l", "--list", dest="list", default=False, action="store_true",
                       help="Prints all formatters available in folder formatters (using duck-typing)")
    parser.add_option ("-s", "--style", dest="styles", default=[], action="append",
                       help="Additional styling options in name:value format (for example separator:\n default is os separator)")
    parser.set_usage ("""%prog [options]""")
    return parser


def parse_args (parser):
    """Parses argument using given parses and check resulting value combination"""
    (options, args) = parser.parse_args ()

    if options.list == True:
        return (options, args)

    if options.input is None:
        print "Error: No input file specified!"
        parser.print_help ()
        sys.exit (1)

    if options.formatter is None:
        print "Error: No formatter specified!"
        parser.print_help ()
        sys.exit (1)

    return (options, args)


def main ():
    parser = create_parser ()
    (options, args) = parse_args (parser)

    # list formatters
    if options.list == True:
        formatters = ProfilerFormatter.list_formatters ()
        print 'Formatter available: \n\t{:s}'.format ('\n\t'.join (formatters))
        sys.exit (0)

    # call main method
    formatter = ProfilerFormatter ()
    result = formatter.convert (options.input, options.output, options.formatter, options.styles)

    # process result
    if result == True:
        sys.exit (0)
    else:
        print result
        sys.exit (1)


# only if this file is main python file, read args
if __name__ == "__main__":
    main ()