#!/usr/bin/python3
# encoding: utf-8
# author:   Jan Brezina, Tomas Krizek

"""
This script check for older version of the IST output format and convert
it into the new format.

Performed actions:
- substitute "type_name" -> "name" (issue of old Record)
- "AbstractRecord" -> "Abstract"

"""

import json
import re
import argparse


def main():
    parser = argparse.ArgumentParser(description='IST normalize script')
    parser.add_argument('input')
    parser.add_argument('output')
    args = parser.parse_args()

    with open(args.input) as f:
        contents = f.read()

    contents = contents.replace('AbstractRecord', 'Abstract')
    contents = contents.replace('type_name', 'name')
    contents = contents.replace('type_full_name', 'full_name')

    data = json.loads(contents)
    if type(data) == list:
        contents = '{\n"ist_nodes" : \n' + contents + '}\n'
    
    with open(args.output, "w") as f:
        f.write(contents)


if __name__ == '__main__':
    main()
