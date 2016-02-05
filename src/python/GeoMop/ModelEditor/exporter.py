"""
usage: exporter.py [-h] [--destination-file [DESTINATION_FILE]] --yaml_file
                   YAML_FILE

Parameters::

    -h, --help
        show this help message and exit
    --destination-file [DESTINATION_FILE]
        The destination file if is different from source file
    --yaml_file YAML_FILE
        YAML input file


Description:
    Script exports new-style YAML Flow123d configuration files into the old-style
    CON files. Constrary to import, no transformations are performed. Export purely
    translates the YAML text into the old CON file.
"""

import os
import sys
__lib_dir__ = os.path.join(os.path.split(
    os.path.dirname(os.path.realpath(__file__)))[0], "common")
sys.path.insert(1, __lib_dir__)

import argparse

from meconfig import cfg


if __name__ == "__main__":
    def main():
        """Launches the export cli."""
        parser = argparse.ArgumentParser(
            description='Export the YAML configuration file to CON format')
        parser.add_argument('--destination-file', nargs='?', default=None,
                            help='The destination file if is different from source file')
        parser.add_argument('--yaml_file', help='YAML input file', required=True)
        args = parser.parse_args()

        if args.destination_file:
            file = args.destination_file
        else:
            file = os.path.splitext(args.yaml_file)[0] + '.con'  # replace extension
        if os.path.isfile(file):
            raise Exception("File already exists")

        cfg.init(None)
        if cfg.open_file(args.yaml_file):
            con_text = cfg.export_file()
            file_d = open(file, 'w')
            file_d.write(con_text)
            file_d.close()

    main()
