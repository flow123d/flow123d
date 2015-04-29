__author__ = 'Jan Hybs'

import pathfix
pathfix.append_to_path ()

import system.versions
system.versions.require_version_2 ()

import sys
from optparse import OptionParser
from profiler.profiler_formatter_module import ProfilerFormatter



def create_parser ():
    parser = OptionParser()
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
    (options, args) = parser.parse_args ()

    if options.list == True:
        return (options, args)

    if options.input == None:
        print "Error: No input file specified!"
        parser.print_help()
        sys.exit(1)

    if options.formatter == None:
        print "Error: No formatter specified!"
        parser.print_help()
        sys.exit(1)

    return (options, args)




def main ():
    parser = create_parser ()
    (options, args) = parse_args (parser)

    # list formatters
    if options.list == True:
        formatters = ProfilerFormatter.list_formatters()
        print 'Formatter available: \n\t{:s}'.format ('\n\t'.join(formatters))
        sys.exit(0)

    # call main method
    formatter = ProfilerFormatter()
    result = formatter.convert (options.input, options.output, options.formatter, options.styles)

    # process result
    if result == True:
        sys.exit(0)
    else:
        print result
        sys.exit(1)


# only if this file is main python file, read args
if __name__ == "__main__" :
    main()