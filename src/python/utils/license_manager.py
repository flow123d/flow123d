# encoding: utf-8
# author:   Jan Hybs

import os
import datetime
from subprocess import CalledProcessError
import re
from utils.logger import Logger


logger = Logger(__name__)


class DoxySection(object):
    def __init__(self, lines=[], formatted=True, section='', obligatory=False):
        self.lines = lines[:]
        self.formatted = formatted
        self.section = section
        self.obligatory = obligatory

    def add(self, line):
        l = line
        if self.formatted:
            l = l.lstrip('* /')
            l = l.rstrip(' \t\r')
        else:
            l = re.sub(r'^\s+\*/?\s?', '', l)

        self.lines.append(l)

    def _fix_lines(self):
        self.lines.reverse()
        for line in self.lines:
            if not line:
                self.lines.remove(line)
            else:
                break
        self.lines.reverse()

    def value(self):
        self._fix_lines()
        lines = self.lines[:]

        if self.formatted:
            if not self.lines and self.obligatory:
                lines = [''] if not self.lines else self.lines
            elif not self.lines and not self.obligatory:
                return ''


        if self.formatted:
            result = '\n * @{:7s} {:s}'.format(self.section, *lines)
            for line in lines[1:]:
                result += '\n *  {space:7s} {line:s}'.format(space='', line=line)
        else:
            result = '\n * @{:7s} '.format(self.section)
            for line in lines:
                result += '\n * {line:s}'.format(line=line)
        return result


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

    def _find_section(self, line):
        l = line
        l = l.lstrip('* /')
        l = l.rstrip(' \t\r')
        simple_section = re.findall(r'@([a-zA-Z]+)[ \t]+(.+)', l)

        if simple_section:
            return simple_section[0][0], simple_section[0][1]

        complex_section = re.findall(r'@([a-zA-Z]+)', l)
        if complex_section:
            return complex_section[0], None
        return None, None

    def add_old_variables(self, variables, old_license):
        # basic variables used in old license
        basic_vars = {
            'brief': DoxySection(section='brief', obligatory=True),
            'ingroup': DoxySection(section='ingroup'),
            'author': DoxySection(section='author'),
            'date': DoxySection(section='date'),
            'section': DoxySection(section='section'),
            'todo': DoxySection(section='todo'),
            'file': DoxySection([os.path.basename(variables['filepath'])], section='file')
        }

        # find old variables
        section = None
        old_vars = { }
        for line in old_license.split('\n'):
            new_section_name, new_section_value = self._find_section(line)
            if new_section_name:
                doxy = DoxySection()
                doxy.section = new_section_name
                doxy.formatted = bool(new_section_value)
                old_vars[new_section_name] = doxy
                section = new_section_name

                if new_section_value:
                    doxy.add(new_section_value)
            elif section:
                old_vars[section].add(line)

        new_vars = basic_vars.copy()
        new_vars.update(old_vars)
        new_vars['brief'].obligatory = True
        keys = new_vars.keys()

        # create _name_ variables
        for key, doxy in new_vars.items():
            new_vars['_' + key + '_'] = doxy.value()
            new_vars[key] = ' '.join(doxy.lines)

        # add them to variable reference
        variables.update(new_vars)


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
