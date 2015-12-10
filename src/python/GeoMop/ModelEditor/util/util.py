"""
Utility classes.

.. codeauthor:: Tomas Krizek <tomas.krizek1@tul.cz>
"""


class TextValue:
    """Represents a value in the input text."""
    def __init__(self, value=None):
        self.value = value
        """the value from input text"""
        self.span = None
        """:class:`.Span` specifies the position of value in input text"""
