"""Workspace selector widget.

..codeauthor:: Tomas Krizek <tomas.krizek1@tul.cz>
"""
import os
import sys

from PyQt5 import QtWidgets
from PyQt5 import QtCore

from .selector_widget import SelectorWidget


class WorkspaceSelectorWidget(SelectorWidget):
    """Widget for selecting the workspace."""

    selected = QtCore.pyqtSignal(str)
    """Signal is emitted when a worksace is selected.

    :param str path: path to the workspace folder
    """

    def __init__(self, parent=None, workspace=None):
        super(WorkspaceSelectorWidget, self).__init__(parent, 'Workspace', workspace)
        self.button.clicked.connect(self.select_workspace)

    @property
    def display_value(self):
        """Convert the workspace variable to a text label."""
        if self._value is None:
            return '(None Selected)'
        # return the folder name
        return os.path.basename(os.path.normpath(self._value))

    def select_workspace(self):
        """Show file dialog to select path to workspace."""
        if self._value is None:
            curr_dir = ''
        else:
            curr_dir = str(os.path.normpath(os.path.join(self._value, '..')))
        sel_dir = QtWidgets.QFileDialog.getExistingDirectory(self, "Choose workspace", curr_dir)
        if not sel_dir:
            sel_dir = None
        elif sys.platform == "win32":
            sel_dir = sel_dir.replace('/', '\\')
        self.value = sel_dir
        self.selected.emit(self._value)


if __name__ == '__main__':
    def main():
        """"Launches widget."""
        import sys
        from PyQt5.QtWidgets import QApplication
        app = QApplication(sys.argv)
        dialog = WorkspaceSelectorWidget(None, '/home/sharp/geo')
        dialog.show()
        dialog.selected.connect(print)
        sys.exit(app.exec_())
    main()

