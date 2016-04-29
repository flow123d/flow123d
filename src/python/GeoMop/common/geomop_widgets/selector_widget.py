"""Selector widget contains a label and a button.

.. codeauthor:: Tomas Krizek <tomas.krizek1@tul.cz>
"""

from PyQt5 import QtWidgets


class SelectorWidget(QtWidgets.QWidget):
    """Widget for selecting a value."""

    def __init__(self, parent=None, label='Selector', value=None):
        super(SelectorWidget, self).__init__(parent)

        self._value = value        
        self.label = QtWidgets.QLabel(label)
        self.button = QtWidgets.QPushButton(self.display_value)
        self.button.setMinimumWidth(180)
        self._layout = QtWidgets.QHBoxLayout()
        self._layout.addWidget(self.label)
        self._layout.addStretch(1)
        self._layout.addWidget(self.button)
        self._layout.setContentsMargins(0,0,0,0)
        
        self.setLayout(self._layout)
        
    @property
    def value(self):
        """Holds the selected value."""
        return self._value

    @value.setter
    def value(self, value):
        """Set the selected value."""
        self._value = value
        self.button.setText(self.display_value)

    @property
    def display_value(self):
        """Holds the value the button displays."""
        return str(self._value)
