# encoding: utf-8
# author:   Jan Hybs

"""
This script is simple iface to license_manager module. Script is used for finding and replacing licenses in given files/locations
Script is expected to run as main python file. Script takes arguments from console and calls methods of module.


Usage: license_manager_script.py [options]

Options:
  -h, --help            show this help message and exit
  -l FILENAME, --license=FILENAME
                        Absolute or relative path to file containing new
                        license, if not specified, license will be only
                        removed
  -f FILENAME, --file=FILENAME
                        Absolute or relative path to file which will be
                        processed
  -d DIRNAME, --dir=DIRNAME
                        Absolute or relative path to dir which will be
                        recursively processed
  -g DIRNAME, --git=DIRNAME
                        Location of git repository (will add additional
                        information to license)
  -s START_SEQ, --start=START_SEQ
                        Start sequence of license, default is '/*!'
  -e END_SEQ, --end=END_SEQ
                        End sequence of license, default is '*/'
  -r, --replace-only    If no license was found, skip file, default is False
  -w, --whitespace      If set, whitespace around license as well as at the
                        beginning of each file will be stripped. After new
                        license two \n will be added
  -o NAME:VALUE, --option=NAME:VALUE
                        Additional formatting variables which will be replaced
                        In license test by using {variable_name} placeholder
                        format syntax you can replace these placeholders
                        dynamically. There are also built in variables:
                        filename, filepath, datetime  You can also format
                        placeholders e.g.: {variable:|^10s} will center value
                        in 10 char e.g. call script with -o foo:bar
                        in license {foo:|^11s}
                        will produce ||||bar||||
                        You can use python format syntax for string only
                        Also date syntax with field datetime e.g.:
                        {datetime:%d:%m:%Y}
                        Additionally if -g or --git is specified other built
                        names will be accessible:
                        last_change, last_author, last_email for each file
                        'branch' globally


## Script license_manager_script

This script is just interface to license_manager module.
Script is used for finding and replacing licenses in given files/locations whilest supporting placeholders.
*For usage run python license_manager_script.py --help*

### Adding locations to script

In order to work user needs to specify files or locations which will be processed.

* Flag -f or --file (can be used multiple times)
  Value should be absolute or relative path to file which will be processed
  
* Flag -d or --dir (can be used multiple times)
  Value should be absolute or relative path to dir which will be recursively processed.
  
### What files will be processed?
All recursively search dirs will process only files with extensions: '.hh', '.cc', '.h', '.c', '.cpp', '.hpp'.
Exception is flag -f or --file which will process given file always. 

### Finding license
User can specify which old license looked like: how it stared and how it ended.
Old license will be regonized as license if it starts at the beginning of the file (whitespace ignored).

* Flag -s or --start (default is /*!)
  Value should be start sequence of the old license
  
* Flag -e or --end (default is */)
  Value should be end sequence of the old license.
  

### New license
User can specify new license which will replace old one. If no license is specified old license will be removed.

* Flag -l or --license default None
  Absolute or relative path to file containing new license
  
### Placeholders
User can specify in license files placeholders which will be processed and replaced. Some placeholders are present by default:

* filepath - absolute path currently processed file (String)
* filename - file name of currently processed file (String)
* datetime - current system date (Date)

if user specifies directory, where is git root, more default values will be available

* Flag -g or --git git root directory:

Additional placeholders available:

* branch - current branch name (String)
* last_change - date of last change (String)
* last_author - name of the last author to make changes (String)
* last_email - email of the last author to make changes (String)

#### Custom placeholders

* Flag -o or --option (can be used multiple times)

  Value should be in form of NAME:VALUE

#### Using placeholders in license file
Placeholders need to be wrapped in braces to be recognized. E.g. {filename}.
Placeholders can be also formated using [this specification](https://docs.python.org/2/library/string.html#format-specification-mini-language).

##### Example 1

file new_license.txt contains following:

/*!
 *
 * Copyright (C) 2007-{datetime:%Y} Technical University of Liberec.  All rights reserved.
 * 
 * @file {filename:>20s} "{filepath}"
 */
 
 
 When called with
 bash
 python license_manager_script.py -f "./path/to/file/main.cc" -l new_license.txt
 
 
 will generate something like this:
 
 /*!
 *
 * Copyright (C) 2007-2015 Technical University of Liberec.  All rights reserved.
 * 
 * @file              main.cc "/usr/jan-hybs/flow123d/path/to/file/main.cc"
 */
 
 ... file content goes here
 
 
 
##### Example 2

file new_license.txt contains following:


/*!
 *
 * Copyright (C) 2007-{datetime:%Y} Technical University of Liberec.  All rights reserved.
 * 
 * @author '{last_author:~>40s} ({last_email})'
 * @file   '{filename:~^40s}'
 * @branch '{branch:~<40s}'
 */
 
 
 When called with
 bash
 python license_manager_script.py -f "./path/to/git/src/main.cc" -l new_license.txt -g "./path/to/git/"
 
 
 will generate something like this:
 
/*!
 *
 * Copyright (C) 2007-2015 Technical University of Liberec.  All rights reserved.
 * 
 * @author '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Jan Hybs (jan.hybs@tul.cz)'
 * @file   '~~~~~~~~~~~~~~~~main.cc~~~~~~~~~~~~~~~~~'
 * @branch 'JHy_ist_formatter~~~~~~~~~~~~~~~~~~~~~~~'
 */
 
 ... file content goes here
 

"""

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

    parser.add_option("-r", "--replace-only", dest="replace_only", default=False, action="store_true",
                      help="If no license was found, skip file, default is %default")

    parser.add_option("-v", "--old-variables", dest="old_variables", default=True, action="store_false",
                      help="Process doxygen variables, default is %default. One-liner variables are in {_name_} format."
                           "{_name_} variables will contain whole formated line.")

    parser.add_option("-w", "--whitespace", dest="whitespace", default=False, action="store_true",
                      help="If set, whitespace around license as well as at the beginning of each "
                           "file will be stripped. After new license two \\n will be added")

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
        variables=variables,
        replace_only=options.replace_only,
        whitespace=options.whitespace,
        old_variables=options.old_variables
    )
    manager.add_locations(options.files, options.dirs)
    manager.add_git(options.git)
    manager.replace_license()


if __name__ == '__main__':
    main()