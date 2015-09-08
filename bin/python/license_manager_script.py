# encoding: utf-8
# author:   Jan Hybs


import sys
import pathfix

pathfix.append_to_path()

import system.versions

system.versions.require_version_2()

from optparse import OptionParser
from utils.license_manager import LicenseManager


def create_parser():
    """Creates command line parse"""
    parser = OptionParser()
    parser.add_option("-l", "--license", dest="license_path", metavar="FILENAME", default=None,
                      help="Absolute or relative path to file containing new license,"
                           " if not specified, license will be only removed")
    parser.add_option("-f", "--file", dest="files", metavar="FILENAME", default=[], action="append",
                      help="Absolute or relative path to file which will be processed")
    parser.add_option("-d", "--dir", dest="dirs", metavar="DIRNAME", default=[], action="append",
                      help="Absolute or relative path to dir which will be recursively processed")

    parser.add_option("-g", "--git", dest="git", metavar="DIRNAME", default=None,
                      help="Location of git repository (will add additional information to license)")

    parser.add_option("-s", "--start", dest="license_start", metavar="START_SEQ", default='/*!',
                      help="Start sequence of license, default is '%default'")

    parser.add_option("-e", "--end", dest="license_end", metavar="END_SEQ", default='*/',
                      help="End sequence of license, default is '%default'")

    parser.add_option("-o", "--option", dest="variables", metavar="NAME:VALUE", default=[], action="append",
                      help="""Additional formatting variables which will be replaced
In license test by using {variable_name} placeholder format syntax
you can replace these placeholders dynamically. There are also built in variables:
filename, filepath, datetime

You can also format placeholders e.g.:
{variable:|^10s} will center value in 10 char
e.g. call script with -o foo:bar
\t\t\t\t\t\t\t\t\t
in license {foo:|^11s}
\t\t\t\t\t\t\t\t\t
will produce ||||bar||||
\t\t\t\t\t\t\t\t\t
You can use python format syntax for string only
\t\t\t\t\t\t\t\t\t
Also date syntax with field datetime e.g.: {datetime:%d:%m:%Y}
\t\t\t\t\t\t\t\t\t
Additionally if -g or --git is specified other built names will be accessible:
\t\t\t\t\t\t\t\t\t
last_change, last_author, last_email for each file
\t\t\t\t\t\t\t\t\t
'branch' globally
                      """)
    parser.set_usage("""%prog [options]""")
    return parser


def parse_args(parser):
    """Parses argument using given parses and check resulting value combination"""
    (options, args) = parser.parse_args()

    if options.license_path is None:
        print "Warning: No license file specified! License will be only removed not replaced!"

    if not options.dirs and not options.files:
        print "Error: No files or dirs specified, nothing to process!"
        parser.print_help()
        sys.exit(1)

    return options, args


def main():
    parser = create_parser()
    (options, args) = parse_args(parser)

    variables = { }
    for name_value in options.variables:
        name, value = name_value.split(':', 1)
        variables[name] = value

    if options.license_path:
        with open(options.license_path, 'r') as fp:
            license_text = fp.read()
    else:
        license_text = ''

    manager = LicenseManager(
        license_text=license_text,
        license_start=options.license_start,
        license_end=options.license_end,
        variables=variables
    )
    manager.add_locations(options.files, options.dirs)
    manager.add_git(options.git)
    manager.replace_license()


if __name__ == '__main__':
    main()