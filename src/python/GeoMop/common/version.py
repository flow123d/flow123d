"""
Module contains version.
"""

import os

__author__ = 'Tomas Krizek'

__git_root_path__ = os.path.join(os.path.split(os.path.dirname(os.path.realpath(__file__)))[0],
                                 '..')

class Version:
    """Geomop version module"""

    def __init__(self):
        """Initializes the class."""

        version_file_path = os.path.join(__git_root_path__, 'VERSION')
        
        try:
            with open(version_file_path) as version_file:
                lines = version_file.readlines()
        except FileNotFoundError:
            lines = []

        self.version = lines[0].strip() if len(lines) > 0 else ''
        """text version geomop pressentation"""
        self.build = lines[1].strip() if len(lines) > 1 else ''
        """text build geomop pressentation"""
        self.date = lines[2].strip() if len(lines) > 2 else ''
        """text date geomop pressentation"""
     
    @staticmethod   
    def get_changelog_path():
        """return changelog path"""
        return os.path.join(__git_root_path__, 'CHANGELOG.md')
