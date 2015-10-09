# encoding: utf-8
# author:   Jan Hybs

import os
import datetime
from subprocess import CalledProcessError
import re
from utils.logger import Logger


logger = Logger(__name__)


class LicenseManager(object):
    """
    Class for adding/removing or replacing license in given files and folders
    """

    def __init__(self, license_text, license_start, license_end,
                 variables={ }, replace_only=False, whitespace=False,
                 old_variables=True):
        self.license_start = license_start
        self.license_stop = license_end
        self.license_text = license_text
        self.variables = variables
        self.replace_only = replace_only
        self.old_variables = old_variables
        self.whitespace = whitespace
        self.extensions = ['.hh', '.cc', '.h', '.c', '.cpp', '.hpp']
        self.files = []
        self.dirs = []
        self.git = None

    def _secure_output(self, command, default=''):
        """
        Run command and return its output or default on error
        :param command:
        :param default:
        """
        from subprocess import check_output

        try:
            return check_output(command, shell=True, cwd=self.git).strip()
        except CalledProcessError as e:
            return default

    def add_git(self, location):
        """
        Add git location
        :param location:
        """
        self.git = location
        if self.git:
            self.variables['branch'] = self._secure_output('git rev-parse --abbrev-ref HEAD', 'unknown-branch')

    def find_license(self, file_content):
        """
        Attempt to find license start and end point
        :param file_content:
        :return: location tuple, or 0,0 tuple
        """
        li_start = file_content.find(self.license_start)
        li_end = file_content.find(self.license_stop)

        if li_start != -1 and li_end != -1:
            li_end += len(self.license_stop)

            if li_start != 0:
                return 0, 0

            return li_start, li_end
        return 0, 0

    def add_locations(self, files=[], dirs=[]):
        """
        Add additional location where to search
        :param files:
        :param dirs:
        """
        self.files.extend(files)
        self.dirs.extend(dirs)

    def replace_license(self):
        """
        Run replacement process for every given file and dir
        :return:
        """
        for rootdir in self.dirs:
            for root, subdirs, files in os.walk(rootdir):
                for filename in files:
                    name, extension = os.path.splitext(filename)
                    if extension.lower() in self.extensions:
                        file_path = os.path.abspath(os.path.join(root, filename))
                        self.process_file(file_path)
        for file_path in self.files:
            self.process_file(file_path)

    def add_old_variables(self, variables, old_license):
        # basic variables used in old license
        basic_vars = {
            'brief': [''],
            'ingroup': [],
            'author': [],
            'date': [],
            'file': [os.path.basename(variables['filepath'])]
        }

        # find old variables
        section = None
        old_vars = { }
        for line in old_license.split('\n'):
            line = line.lstrip('* /')
            line = line.rstrip('\r')
            if not line.strip():
                section = None
                continue

            result = re.findall(r'@([a-zA-Z]+)[ \t]+(.+)', line)
            if result:
                var_name, var_value = result[0]
                var_value = var_value.strip()
                old_vars[var_name] = [var_value]
                section = var_name
            elif section:
                old_vars[var_name].append(line)

        # max_length = 7#len(max(old_vars.keys(), key=len)) + 1
        new_vars = basic_vars.copy()
        new_vars.update(old_vars)
        keys = new_vars.keys()

        # create _name_ variables
        for key, values in new_vars.items():
            if not values:
                new_vars['_' + key + '_'] = ''
                continue

            if values[0].strip() == '???':
                new_vars['_' + key + '_'] = ' '
                continue

            fmt = '\n * @{:7s} {:s}' + '\n *  {space:7s} {:s}'.join([''] * len(values))
            new_vars['_' + key + '_'] = fmt.format(key, *values, space='')

        # add them to variable reference
        variables.update(new_vars)

        # convert lists to str
        for key in keys:
            variables[key] = '\n'.join(variables[key])

    def process_file(self, file_path):
        """
        Process single file
        :param file_path:
        :return:
        """
        logger.info ('Processing {:s}'.format(file_path))
        with open(file_path, 'r') as fp:
            file_content = fp.read()

        file_content = file_content.lstrip()
        li_start, li_end = self.find_license(file_content)

        if li_start == 0 and li_end == 0:
            if self.replace_only:
                logger.info('File {:s} skipped, no license found'.format(file_path))
                return

        if li_start is not None and li_end is not None:
            old_license = file_content[li_start:li_end]
            file_content = file_content[li_end:].lstrip()
            variables = self.variables.copy()
            variables.update(
                {
                    'filepath': file_path,
                    'filename': os.path.basename(file_path),
                    'datetime': datetime.datetime.now()
                })

            if self.old_variables:
                self.add_old_variables(variables, old_license)

            if self.git:
                variables.update({
                    'last_change': self._secure_output('git log -1 --format=%cd '+ file_path),
                    'last_author': self._secure_output('git log -1 --format=%cn '+ file_path),
                    'last_email':  self._secure_output('git log -1 --format=%ce '+ file_path)
                })

            license_text = self.license_text.strip() if self.whitespace else self.license_text
            license_text = license_text.format(**variables)

            with open(file_path, 'w') as fp:
                if license_text:
                    fp.write(license_text)
                    if self.whitespace:
                        fp.write('\n\n')
                fp.write(file_content)
