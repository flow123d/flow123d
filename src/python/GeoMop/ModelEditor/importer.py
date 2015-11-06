"""
usage:
    importer.py [-h]
    [--transformation-name [TRANSFORMATION_NAME [TRANSFORMATION_NAME ...]]]
    [--destination-file [DESTINATION_FILE]] --con_file CON_FILE

Parameters::

    -h, --help
        show this help message and exit
    --transformation-name [TRANSFORMATION_NAME [TRANSFORMATION_NAME ...]]
        Transformation rules contained in the Model Editor that are processed after import
    --destination-file [DESTINATION_FILE]
        The destination file if is different from source file
    --con_file CON_FILE
        Con input file

Description:
    Importer translate Flow123d con file to yaml format. TYPE parameters is
    transformed to tags (!value). Importer try rewrite references to yaml
    format and place it to appropriate place. Comments are copied and
    place to yaml file.  If is set some transformation rules, they are processed
    after import in set order by :mod:`data.yaml.transformator`.
"""

import os
import sys
__lib_dir__ = os.path.join(os.path.split(
    os.path.dirname(os.path.realpath(__file__)))[0], "common")
sys.path.insert(1, __lib_dir__)

from data.meconfig import MEConfig as cfg
import argparse

if __name__ == "__main__":
    def main():
        """Launches the import cli."""
        parser = argparse.ArgumentParser(
            description='Import the yaml configarion file from con format')
        parser.add_argument('--transformation-name', nargs='*',
                            help='Transformation rules contained in the Model Editor that are '
                                 'processed after import')
        parser.add_argument('--destination-file', nargs='?', default=None,
                            help='The destination file if is different from source file')
        parser.add_argument('--con_file', help='Con input file', required=True)
        args = parser.parse_args()

        if args.destination_file:
            file = args.destination_file
        else:
            if args.con_file[:-4]:
                file = args.con_file[:-4] + ".yaml"
            else:
                file = args.con_file[:-4] + ".yaml"
        if os.path.isfile(file):
            raise Exception("File already exists")

        cfg.init(None)
        cfg.import_file(args.con_file)
        if args.transformation_name is not None:
            for transf in args.transformation_name:
                cfg.transform(transf)

        file_d = open(file, 'w')
        file_d.write(cfg.document)
        file_d.close()

    main()
