"""
Library for work with state of GeoMops applications
"""

import os
import yaml

if 'APPDATA' in os.environ:
    __config_dir__ = os.path.join(os.environ['APPDATA'], 'GeoMop')
else:
    __config_dir__ = os.path.join(os.environ['HOME'], '.geomop')

try:
    if not os.path.isdir(__config_dir__):
        os.makedirs(__config_dir__)
except:
    raise Exception('Cannot create config directory')


def get_config_file(name, directory=None):
    """
    Get config object from file Name.cfg in config directory

    return: Config object or None (if file not exist)
    """
    if directory is not None:
        directory = os.path.join(__config_dir__, directory)
        if not os.path.isdir(directory):
            return None
    else:
        directory = __config_dir__
    file_name = os.path.join(directory, name+'.yaml')
    try:
        yaml_file = open(file_name, 'r')
    except (FileNotFoundError, IOError):
        return None
    config = yaml.load(yaml_file)
    yaml_file.close()
    return config


def save_config_file(name, config, directory=None):
    """Save config object to file Name.cfg in config directory"""
    if directory is not None:
        directory = os.path.join(__config_dir__, directory)
        try:
            if not os.path.isdir(directory):
                os.makedirs(directory)
        except:
            raise Exception('Cannot create config directory: ' + directory)
    else:
        directory = __config_dir__
    file_name = os.path.join(directory, name+'.yaml')
    yaml_file = open(file_name, 'w')
    yaml.dump(config, yaml_file)
    yaml_file.close()


def delete_config_file(name, directory=None):
    """
    Delete config file Name.cfg from config directory
     """
    if directory is not None:
        directory = os.path.join(__config_dir__, directory)
        if not os.path.isdir(directory):
            return
    else:
        directory = __config_dir__
    file_name = os.path.join(directory, name+'.yaml')
    try:
        os.remove(file_name)
    except  (RuntimeError, IOError):
        return
