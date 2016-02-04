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
import json
import pathfix
pathfix.append_to_path()

import system.versions
system.versions.require_version_2()

import sys
from optparse import OptionParser
from utils.logger import Logger

from ist.nodes import TypeRecord, TypeAbstract, TypeSelection, TypeString, TypeDouble, TypeInteger, TypeBool, TypeArray, \
    TypeParameter, TypeFilename


# all acceptable input_type values
registered_nodes = {
    'Record': TypeRecord,
    'AbstractRecord': TypeAbstract,
    'Abstract': TypeAbstract,
    'Selection': TypeSelection,
    'String': TypeString,
    'Double': TypeDouble,
    'Integer': TypeInteger,
    'FileName': TypeFilename,
    'Bool': TypeBool,
    'Array': TypeArray,
    'Parameter': TypeParameter
}


def create_parser():
    """Creates command line parse"""
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="input", metavar="FILENAME", default=None,
                      help="Absolute or relative path to JSON file which will be processed")
    parser.add_option("-o", "--output", dest="output", metavar="FILENAME", default=None,
                      help="Absolute or relative path output file which will be generated/overwritten")
    parser.add_option("-f", "--format", dest="format", metavar="FORMAT", default="tex",
                      help="tex|html output format")
    parser.add_option("--log", dest="loglevel", metavar="LEVEL", default="warning",
                      help="Logging level default is %default, options are (debug, info, warning, error, critical)")
    parser.set_usage("""%prog [options]""")
    return parser


def parse_args(parser):
    """Parses argument using given parses and check resulting value combination"""
    (options, args) = parser.parse_args()

    if options.input is None:
        Logger.instance().warning("Error: No input file specified!")
        parser.print_help()
        sys.exit(1)

    if options.output is None:
        Logger.instance().warning("Error: No output file specified!")
        parser.print_help()
        sys.exit(1)

    return options, args


def main():
    parser = create_parser()
    options, args = parse_args(parser)

    # create instance of formatter
    from ist.ist_formatter_module import ISTFormatter
    formatter = ISTFormatter()

    # read input json file
    with file(options.input, 'r') as fp:
        json_data = json.load(fp)
        json_data = json_data['ist_nodes'] if 'ist_nodes' in json_data else json_data

        # filter out unsupported types, they won't be formatted
        items = list()
        for json_item in json_data:
            input_type = json_item['input_type'] if 'input_type' in json_item else None
            if input_type in registered_nodes:

                item = registered_nodes[input_type]()
                item.parse(json_item)
                items.append(item)
            else:
                Logger.instance().info(' - item type not supported: %s' % str(item))

    # convert to tex format
    if options.format.lower() in ('tex', 'latex'):
        Logger.instance().info('Formatting ist to tex format')
        formatter.json2latex(items, options.output)
        sys.exit(0)

    # convert to HTML format
    if options.format.lower() in ('html', 'html5', 'www', 'htm'):
        Logger.instance().info('Formatting ist to html format')
        formatter.json2html(items, options.output)
        sys.exit(0)

    Logger.instance().error("Error: Unsupported format '{:s}'".format(options.format))
    sys.exit(1)


# only if this file is main python file, read args
if __name__ == "__main__":
    main()