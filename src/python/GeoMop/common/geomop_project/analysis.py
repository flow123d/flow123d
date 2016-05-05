"""Analysis - class + load/save functions.

.. codeauthor:: Tomas Krizek <tomas.krizek1@tul.cz>
"""

import json

ANALYSIS_FILE_EXT = 'data'
ANALYSIS_DEFAULT_NAME = 'analysis'


class Analysis:

    def __init__(self, name=ANALYSIS_DEFAULT_NAME, files=None, params=None, **kwargs):
        self.name = name
        """name of the analysis"""
        self.files = []
        """a list of files used in analysis"""
        self.params = {}
        """key: value pairs of parameter values"""
        if files is not None and isinstance(files, list):
            self.files = files
        if params is not None and isinstance(params, dict):
            self.params = params

    @property
    def filename(self):
        """filename of the analysis"""
        return self.name + '.' + ANALYSIS_FILE_EXT


def load_analysis_file(file_path):
    """Load the analysis files."""
    with open(file_path) as file:
        data = json.load(file)
        analysis = Analysis(**data)
    return analysis


def save_analysis_file(file_path, analysis):
    """Load the analysis files."""
    with open(file_path, 'w') as file:
        json.dump(analysis.__dict__, file)
