"""
.. codeauthor:: Tomas Krizek <tomas.krizek1@tul.cz>
"""

import os

from geomop_util import Serializable
import config

from .analysis import (ANALYSIS_FILE_EXT, save_analysis_file, load_analysis_file,
                       ANALYSIS_DEFAULT_NAME)


PROJECT_MAIN_FILE_EXT = 'project'
PROJECT_MAIN_FILE = '.'+PROJECT_MAIN_FILE_EXT


class InvalidProject(Exception):
    pass


class Parameter:
    """A parameter in a config file."""
    def __init__(self, name, type=None):
        self.name = name
        self.type = type


class File:
    """Represents a file entry in a config file."""
    def __init__(self, file_path, params=None):
        self.file_path = file_path
        if params is None:
            self.params = []
        else:
            self.params = params


class Project:
    """Project settings and data."""

    __serializable__ = Serializable(
        excluded=['filename',
                  'workspace',
                  'name',
                  'project_dir'],
        composite={'params': Parameter,
                   'files': File}
    )

    current = None
    """currently opened project"""

    def __init__(self, filename=None, **kwargs):
        self.filename = filename
        self.workspace = None
        self.name = None
        self.params = kwargs['params'] if 'params' in kwargs else []
        self.files = kwargs['files'] if 'files' in kwargs else []
        self.project_dir = ''

    @staticmethod
    def open(workspace, project_name):
        """Retrieve project from settings by its name and workspace."""
        directory = os.path.join(workspace, project_name)
        project = config.get_config_file('',
                                         directory=directory,
                                         extension=PROJECT_MAIN_FILE_EXT,
                                         cls=Project)
        if project is None:
            raise InvalidProject("Current project is invalid.")

        project.project_dir = directory
        project.filename = os.path.join(directory, PROJECT_MAIN_FILE)
        project.workspace = workspace
        project.name = project_name
        return project

    @staticmethod
    def exists(workspace, project_name):
        """Determine whether project exists in a workspace."""
        if not workspace or not project_name:
            return False
        project_filename = os.path.join(workspace, project_name, PROJECT_MAIN_FILE)
        if not os.path.isfile(project_filename):
            return False
        return True

    @staticmethod
    def notify(data):
        """Observer method to update current project."""
        if data.workspace is None or data.project is None:
            Project.current = None
            return
        if (Project.current is None or data.workspace != Project.current.workspace or
                data.project != Project.current.name):
            if Project.exists(data.workspace, data.project):
                project = Project.open(data.workspace, data.project)
                Project.current = project
            else:
                data.project = None
                Project.current = None

    @staticmethod
    def reload_current():
        """Read the current project file and updated the project data."""
        if Project.current is None:
            return None
        project = Project.open(Project.current.workspace, Project.current.name)
        Project.current = project
        return project

    def save(self):
        """Save the current project as a config file."""
        config.save_config_file('', self, self.project_dir, PROJECT_MAIN_FILE_EXT)

    def load_analysis(self, name=ANALYSIS_DEFAULT_NAME):
        """Loads a project analysis by its name."""
        file_path = os.path.join(self.project_dir, name + '.' + ANALYSIS_FILE_EXT)
        try:
            return load_analysis_file(file_path)
        except Exception:
            return None

    def save_analysis(self, analysis):
        """Save the analysis within the current project."""
        file_path = os.path.join(self.project_dir, analysis.name + '.' + ANALYSIS_FILE_EXT)
        save_analysis_file(file_path, analysis)

    def get_all_analyses(self):
        """Return all project analyses."""
        analyses = []
        for file_name in os.listdir(self.project_dir):
            if file_name.endswith('.' + ANALYSIS_FILE_EXT):
                analysis = self.load_analysis(file_name[:-len('.' + ANALYSIS_FILE_EXT)])
                analyses.append(analysis)
        return analyses

    def get_current_analysis(self):
        """Return current analysis."""
        analyses = self.get_all_analyses()
        for analysis in analyses:
            if analysis.name == ANALYSIS_DEFAULT_NAME:
                return analysis
        else:
            return analyses[0] if len(analyses) > 0 else None

    def make_relative_path(self, file_path):
        """Makes the path relative to project_dir."""
        if not self.project_dir:
            return file_path
        if not file_path.startswith(self.project_dir):
            # assume file_path is already relative to project
            return file_path
        return file_path[len(self.project_dir)+1:]

    def merge_params(self, params):
        """Merge another param collection into this one."""
        for new_param in params:
            # check if param doesn't exist, then add it
            exists = False
            for curr_param in self.params:
                if curr_param.name == new_param.name:
                    exists = True
                    break
            if not exists:
                self.params.append(new_param)

    def is_abs_path_in_project_dir(self, file_path):
        """Whether file exists in the project directory or subdirectories."""
        if not file_path or not self.project_dir:
            return False
        return file_path.startswith(self.project_dir)

    def add_file(self, file_path, params=None):
        """Add or sync a file. If file exists, update parameters."""
        file_path = self.make_relative_path(file_path)
        for file in self.files:
            if file_path == file.file_path:
                if params is not None:
                    file.params = params
                return
        # file is not registered
        file = File(file_path, params)
        self.files.append(file)
        return

