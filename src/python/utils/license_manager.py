# encoding: utf-8
# author:   Jan Hybs

import os
import datetime
from subprocess import CalledProcessError


class LicenseManager(object):
    """
    Class for adding/removing or replacing license in given files and folders
    """
    def __init__(self, license_text, license_start, license_end, variables={}):
        self.license_start = license_start
        self.license_stop = license_end
        self.license_text = license_text
        self.variables = variables
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

    def process_file(self, file_path):
        """
        Process single file
        :param file_path:
        :return:
        """
        print ('Processing {:s}'.format(file_path))
        with open(file_path, 'r') as fp:
            file_content = fp.read()

        file_content = file_content.lstrip()
        li_start, li_end = self.find_license(file_content)

        if li_start is not None and li_end is not None:
            file_content = file_content[li_end:].lstrip()
            variables = self.variables.copy()
            variables.update(
                {
                    'filepath': file_path,
                    'filename': os.path.basename(file_path),
                    'datetime': datetime.datetime.now()
                })

            if self.git:
                variables.update({
                    'last_change': self._secure_output('git log -1 --format=%cd '+ file_path),
                    'last_author': self._secure_output('git log -1 --format=%cn '+ file_path),
                    'last_email':  self._secure_output('git log -1 --format=%ce '+ file_path)
                })

            license_text = self.license_text.strip()
            license_text = license_text.format(**variables)

            with open(file_path, 'w') as fp:
                if license_text:
                    fp.write(license_text)
                    fp.write('\n\n')
                fp.write(file_content)
