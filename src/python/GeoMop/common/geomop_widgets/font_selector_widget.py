"""Font selector widget.

..codeauthor:: Tomas Krizek <tomas.krizek1@tul.cz>
"""
import os

from PyQt5 import QtWidgets
from PyQt5 import QtCore
from PyQt5 import QtGui

from .selector_widget import SelectorWidget


class FontSelectorWidget(SelectorWidget):
    """Widget for selecting the font."""

    def __init__(self, parent=None, font_string=None):
        value = QtGui.QFont()
        value.fromString(font_string)
        super(FontSelectorWidget, self).__init__(parent, 'Font', value)
        self.button.clicked.connect(self.select_font)

    @property
    def display_value(self):
        """Convert the font variable to a text label."""
        if self._value is None:
            return '(None Selected)'
        # return the font name and size
        return ','.join(self._value.toString().split(',')[:2])  # don't display extra params

    def select_font(self):
        """Show font select dialog."""
        selected_font, ok = QtWidgets.QFontDialog.getFont(
            self._value, self, options=QtWidgets.QFontDialog.MonospacedFonts)
        if ok:
            self.value = selected_font
