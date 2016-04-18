#!/usr/bin/python
# -*- coding: utf-8 -*-
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

from __future__ import absolute_import

import pathfix
pathfix.append_to_path()

import system.versions
system.versions.require_version_2()

import os, sys, json
from optparse import OptionParser
from utils.logger import Logger
from ist.base import InputType
from ist.utils.htmltree import htmltree
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
        ist_info = {
            'version': json_data['version']['flow123d_version'] if 'version' in json_data else 'Input reference'
        }
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

    # if we have all items parsed we create references
    for item in items:
        if getattr(item, 'input_type', InputType.UNKNOWN) == InputType.MAIN_TYPE:
            if item.input_type == InputType.RECORD:
                for key in getattr(item, 'keys', []):
                    if key.type.get_reference().input_type == InputType.ARRAY:
                        key.type.get_reference().subtype.get_reference().add_ref(item)
                    else:
                        key.type.get_reference().add_ref(item)
            if item.input_type == InputType.ABSTRACT_RECORD:
                for imp in getattr(item, 'implementations', []):
                    imp.get_reference().add_ref(item)

    # sort items by type and name
    items = sorted(items, key=lambda x: '{}{}'.format(x.input_type.value, x.name))

    # convert to tex format
    if options.format.lower() in ('tex', 'latex'):
        Logger.instance().info('Formatting ist to tex format')
        formatter.json2latex(items, options.output, info=ist_info)
        if os.path.isfile(options.output):
            print 'Ok: File "{:s}" created'.format(options.output)
            sys.exit(0)
        else:
            print 'Error: File "{:s}" does not exists'.format(options.output)
            sys.exit(1)

    # convert to HTML format
    if options.format.lower() in ('html', 'html5', 'www', 'htm'):
        Logger.instance().info('Formatting ist to html format')
        formatter.json2html(items, options.output, info=ist_info)
        if os.path.isfile(options.output):
            print 'Ok: File "{:s}" created'.format(options.output)
            sys.exit(0)
        else:
            print 'Error: File "{:s}" does not exists'.format(options.output)
            sys.exit(1)

    if options.format.lower() in ('markdown', 'md'):
        Logger.instance().info('Testing markdown')
        text = '''
# Using markdown in description

**Description field** supports markdown syntax (support is partial and some techniques may not work in Python markdown implementation).

## Links

Link to record [[root]] selection [[DG_output_fields]] or abstract [[Transport]]. All links are in the same format.
If `link_name` is specified in `attributes` (let say DG_output_fields has link_name of DG), we can use that [[DG]]
Record and Selection types offer links to their keys/values, so you can write [[DG#porosity]], to specify link text use following syntax [[DG#porosity:poro]]

or link to key in Root item [[root#flow123d_version]]




Every name should be unique, if `link_name` is duplicate first occurrence will be used!
To avoid conflict with `link_name` or `name` we can use type specification like this [[record#root]]. 3 types are registered:

 1. type RECORD supporting prefixes:

   - r
   - record

 2. type SELECTION supporting prefixes:

   - s
   - selection

 3. type ABSTRACT supporting prefixes:

   - a
   - ar
   - abstract


## Basics
We can write **bold** statements (or *italic* if needed). We can also ~~strkkethrg~~ strikethrough some text to express some change.

Another usage can be in lists, we can write both unordered and ordered list. Important is to place one empty line before list starts.
Unordered list:

 - important item
 - another important item

Ordered list have same rules:

 1. item number 1
 2. and item number 2

To write code section with monospaced font, such as variable use `this` syntax.

**Note** Use line breaks \\n chars sparely, and only break text-flow if necessarily.

Full markdown specification can be found [here](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet) but since Python is used to parse markdown, there may be slight differences.
(($q_w$))
        '''
        o = htmltree()
        o.description(text)
        print o.dump()
        sys.exit(0)

    Logger.instance().error("Error: Unsupported format '{:s}'".format(options.format))
    sys.exit(1)


# only if this file is main python file, read args
if __name__ == "__main__":
    main()