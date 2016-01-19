"""GeoMop's dialogs shared between different contexts.

Dialogs have a standard appearance across contexts.

.. codeauthor:: Pavel Richter <pavel.richter@tul.cz>
"""
import PyQt5.QtWidgets as QtWidgets


class GMErrorDialog(QtWidgets.QMessageBox):
    """Error dialog for GeoMop graphic applications"""
    def __init__(self, parent):
        """Initialize standard QDialog with GeoMop specific properties."""
        super(GMErrorDialog, self).__init__(parent)
        self.setIcon(QtWidgets.QMessageBox.Critical)

    def open_error_dialog(self, text, error=None, title="GeoMop Error"):
        """Display dialog with title, text and error message in detail."""
        self.setWindowTitle(title)
        self.setText(text)
        if error is None:
            self.setInformativeText("")
        else:
            self.setInformativeText(type(error).__name__ + ": " + str(error))
        super(GMErrorDialog, self).exec_()
        
    def open_error_report_dialog(self, errors, msg="Process notifications:" ,  title="GeoMop Operation Report"):
        """Display dialog with title, text and error message in detail."""
        self.setWindowTitle(title)
        self.setInformativeText("\n".join(errors))
        self.setText(msg)
        self.setIcon(QtWidgets.QMessageBox.Information)
        super(GMErrorDialog, self).exec_()    
