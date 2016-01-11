"""
Module contains version.

.. codeauthor:: Tomas Krizek <tomas.krizek1@tul.cz>
"""

import os

__root_path__ = (os.path.split(os.path.dirname(os.path.realpath(__file__)))[0])


def get_root_file_path(filename):
    """File path to a file in a root folder."

    Git repo has has a different place for root than installed version of app.
    """
    file_path = os.path.join(__root_path__, filename)  # installed on Win
    if not os.path.isfile(file_path):
        file_path = os.path.join(__root_path__, '..', filename)  # for git repo
    return file_path


class Version:
    """Geomop version module"""

    def __init__(self):
        """Initializes the class."""
        try:
            with open(Version.VERSION_FILE_PATH) as version_file:
                lines = version_file.readlines()
        except FileNotFoundError:
            lines = []

        self.version = lines[0].strip() if len(lines) > 0 else ''
        """text version geomop pressentation"""
        self.build = lines[1].strip() if len(lines) > 1 else ''
        """text build geomop pressentation"""
        self.date = lines[2].strip() if len(lines) > 2 else ''
        """text date geomop pressentation"""

    CHANGELOG_FILE_PATH = get_root_file_path('CHANGELOG.md')
    VERSION_FILE_PATH = get_root_file_path('VERSION')

