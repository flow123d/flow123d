"""Dialog for creating new projects.

.. codeauthor:: Tomas Krizek<tomas.krizek1@tul.cz>
"""
import os

from PyQt5 import QtWidgets

from geomop_project import PROJECT_MAIN_FILE


class CreateProjectDialog(QtWidgets.QDialog):
    """Dialog for creating projects within a workspace."""

    def __init__(self, parent, config, title='New Project'):
        """Initialize the class."""
        super(CreateProjectDialog, self).__init__(parent)
        self.config = config
        name_label = QtWidgets.QLabel('Project Name')
        self.name_field = QtWidgets.QLineEdit()
        button_box = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok |
            QtWidgets.QDialogButtonBox.Cancel
        )
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        self.ok_button = button_box.button(QtWidgets.QDialogButtonBox.Ok)
        name_layout = QtWidgets.QHBoxLayout()
        name_layout.addWidget(name_label)
        name_layout.addWidget(self.name_field)
        main_layout = QtWidgets.QVBoxLayout()
        main_layout.addLayout(name_layout)
        main_layout.addWidget(button_box)
        self.setLayout(main_layout)
        self.setWindowTitle(title)

    def accept(self):
        """Handle a confirmation."""
        name = self.name_field.text()
        path = os.path.join(self.config.workspace, name)
        if os.path.exists(path):
            QtWidgets.QMessageBox.critical(
                self, 'Name is not unique',
                "Can't create project. The selected project name already exists."
            )
            return
        os.mkdir(path)
        open(os.path.join(path, PROJECT_MAIN_FILE), 'w').close()
        self.config.project = name
        self.config.save()
        return super(CreateProjectDialog, self).accept()

