"""Module contains an edit menu widget.

.. codeauthor:: Tomas Krizek <tomas.krizek1@tul.cz>
"""
import os

from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QMenu, QAction

from geomop_dialogs import CreateProjectDialog
from geomop_project import PROJECT_MAIN_FILE


class ProjectMenu(QMenu):
    """Menu with projects."""

    def __init__(self, parent, config, title='&Project'):
        """Initializes the class."""
        super(ProjectMenu, self).__init__(parent)
        self.config = config
        self._group = QtWidgets.QActionGroup(self, exclusive=True)
        self._group.triggered.connect(self._project_selected)
        self._actions = []
        self._new_project_action = QAction('&New Project', self)
        self._new_project_action.triggered.connect(self._create_new_project)
        self.addAction(self._new_project_action)
        self.addSeparator()
        self.aboutToShow.connect(self.reload_projects)
        self.setTitle(title)

    def reload_projects(self):
        """Reload all projects from current workspace."""
        # remove old actions
        for action in self._actions:
            self.removeAction(action)
            self._group.removeAction(action)
        self._actions = []
        # find projects in workspace
        if not self.config.workspace:
            return
        for name in os.listdir(self.config.workspace):
            path = os.path.join(self.config.workspace, name)
            if is_project(path):
                action = QtWidgets.QAction(name, self, checkable=True)
                action.setChecked(self.config.project == name)
                self.addAction(action)
                self._group.addAction(action)
                self._actions.append(action)

        for action in self._actions:
            self.addAction(action)

    def _project_selected(self):
        """Handle project selection."""
        action = self._group.checkedAction()
        self.config.project = action.text()
        self.config.save()

    def _create_new_project(self):
        """Handle creation of a new project."""
        if not self.config.workspace:
            QtWidgets.QMessageBox.information(
                self, 'No Workspace',
                'Please select workspace in Settings before creating a project.')
        else:
            CreateProjectDialog(self, self.config).show()


def is_project(path):
    """Determine if a path is a project directory."""
    if os.path.isdir(path):
        for file_ in os.listdir(path):
            if file_.lower() == PROJECT_MAIN_FILE:
                return True
    return False

