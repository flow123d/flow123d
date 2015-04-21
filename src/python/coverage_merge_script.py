__author__ = 'Jan Hybs'


import sys
from optparse import OptionParser
from coverage import coverage_merge_module


# parse arguments
def create_parser ():
    newline = 10*'\t';
    parser = OptionParser(usage="%prog [options] [file1 file2 ... filen]", version="%prog 1.0",
        epilog = "If no files are specified all xml files in current directory will be selected. \n" +
                 "Useful when there is not known precise file name only location")

    parser.add_option("-o", "--output",     dest="filename",    default="coverage-merged.xml",
        help="output file xml name", metavar="FILE")
    parser.add_option("-p", "--path",       dest="path",        default="./",
        help="xml location, default current directory", metavar="FILE")
    parser.add_option("-l", "--log",        dest="loglevel",    default="DEBUG",
        help="Log level DEBUG, INFO, WARNING, ERROR, CRITICAL")
    parser.add_option("-f", "--filteronly", dest="filteronly",  default=False, action='store_true',
        help="If set all files will be filtered by keep rules otherwise "+
             "all given files will be merged and filtered.")
    parser.add_option("-s", "--suffix",     dest="suffix",      default='',
        help="Additional suffix which will be added to filtered files so they original files can be preserved")
    parser.add_option("-k", "--keep",       dest="packagefilters", default=None,  metavar="NAME", action="append",
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
    (options, args) = parser.parse_args ()

    # for now, no check needed

    return (options, args)


if __name__ == '__main__':
    parser = create_parser ()
    (options, args) = parse_args (parser)

    coverage_merge_module.merge (options, args)

