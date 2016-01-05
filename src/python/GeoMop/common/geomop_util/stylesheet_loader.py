"""
Module contains helper functions for loading stylesheets.

.. codeauthor:: Tomas Krizek <tomas.krizek1@tul.cz>
"""

import os
import platform


def load_stylesheet(path, resource_dir):
    """Loads the stylesheet file to a string.

    This loading function fixes the paths in css files. All files in the css are considered
    relative to the resources directory. This function prepends all url(...) with the absolute
    path to resource directory. This makes the css resources load correctly independently of
    the working directory.

    :param str path: absolute path to the stylesheet
    :param str resource_dir: absolute path to resource directory
    :return: stylesheet loaded as text
    :rtype: str
    """
    stylesheet = ''

    if os.path.isfile(path):
        with open(path) as file_:
            stylesheet = file_.read()

        if platform.system() == 'Linux':
            stylesheet = stylesheet.replace("url('resources/", "url('" + resource_dir + os.path.sep)

    return stylesheet
