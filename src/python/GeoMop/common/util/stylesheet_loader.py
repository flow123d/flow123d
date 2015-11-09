"""
Module contains helper functions for loading stylesheets.
"""

import os

__author__ = 'Tomas Krizek'

__stylesheet_dir__ = os.path.join(os.getcwd(), 'resources', 'css') + os.path.sep


def load_stylesheet(name):
    """Loads the stylesheet file to a string."""
    if not name.endswith('.css'):
        name += '.css'

    path = __stylesheet_dir__ + name
    if os.path.isfile(path):
        with open(path) as file_:
            return file_.read()

    return ''
