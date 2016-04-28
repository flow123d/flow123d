"""Library for work with state of GeoMops applications.

.. codeauthor:: Pavel Richter <pavel.richter@tul.cz>
"""

import os
import yaml

from geomop_util import Serializable


if 'APPDATA' in os.environ:
    __config_dir__ = os.path.join(os.environ['APPDATA'], 'GeoMop')
else:
    __config_dir__ = os.path.join(os.environ['HOME'], '.geomop')

try:
    if not os.path.isdir(__config_dir__):
        os.makedirs(__config_dir__)
except:
    raise Exception('Cannot create config directory')


def get_config_file(name, directory=None, cls=None, extension='yaml'):
    """
    Get config object from filename in config directory

    return: Config object or None (if file not exist)
    """
    if directory is not None:
        directory = os.path.join(__config_dir__, directory)
        if not os.path.isdir(directory):
            return None
    else:
        directory = __config_dir__
    file_name = os.path.join(directory, name+'.'+extension)
    try:
        yaml_file = open(file_name, 'r')
    except (FileNotFoundError, IOError):
        return None
    config = yaml.load(yaml_file)
    yaml_file.close()
    config = Serializable.load(config, cls)
    return config


def save_config_file(name, config, directory=None, extension='yaml'):
    """Save config object to file name.extension in config directory"""
    if directory is not None:
        directory = os.path.join(__config_dir__, directory)
        try:
            if not os.path.isdir(directory):
                os.makedirs(directory)
        except:
            raise Exception('Cannot create config directory: ' + directory)
    else:
        directory = __config_dir__
    file_name = os.path.join(directory, name+'.'+extension)
    data = Serializable.dump(config)
    yaml_file = open(file_name, 'w')
    yaml.dump(data, yaml_file)
    yaml_file.close()


def delete_config_file(name, directory=None, extension='yaml'):
    """
    Delete config file name.extension from config directory
     """
    if directory is not None:
        directory = os.path.join(__config_dir__, directory)
        if not os.path.isdir(directory):
            return
    else:
        directory = __config_dir__
    file_name = os.path.join(directory, name+'.'+extension)
    try:
        os.remove(file_name)
    except (RuntimeError, IOError):
        return
