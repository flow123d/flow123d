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

import argparse

from meconfig import cfg


def do_transform(con_file, yaml_file, dest_file, transformation_name):
        cfg.init(None)        
        if con_file:
            cfg.import_file(con_file)
            source_file = con_file
        elif yaml_file:
            cfg.open_file(yaml_file)
            source_file = yaml_file

        # check destination file
        if dest_file:
            file = dest_file
        else:
            file = os.path.splitext(source_file)[0] + '.yaml'  # replace extension
        if os.path.isfile(file):
            raise Exception("File already exists")

        # apply transformations        
        if transformation_name is not None:
            for transf in transformation_name:                
                cfg.transform(transf)

        file_d = open(file, 'w')
        file_d.write(cfg.document)
        file_d.close()
    
    


if __name__ == "__main__":
    def main():
        """Launches the import cli."""
        parser = argparse.ArgumentParser(
            description='Import the YAML configuration file from CON format')
        parser.add_argument('--transformation-name', nargs='*',
                            help='Transformation rules contained in the Model Editor that are '
                                 'processed after import')
        parser.add_argument('--destination-file', nargs='?', default=None,
                            help='The destination file if is different from source file')
        parser.add_argument('--con_file', help='CON input file', nargs='?')
        parser.add_argument('--yaml_file', help='YAML input file', nargs='?')
        args = parser.parse_args()

        do_transform(args.con_file, args.yaml_file, args.destination_file, args.transformation_name)

    # call main
    main()
