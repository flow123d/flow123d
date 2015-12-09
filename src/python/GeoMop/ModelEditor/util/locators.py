"""
Classes for determining node and cursor location.

.. codeauthor:: Tomas Krizek <tomas.krizek1@tul.cz>
"""

from geomop_util import ComparableMixin


class Position(ComparableMixin):
    """Marks a cursor position in text."""
    def __init__(self, line=None, column=None):
        self.line = line
        """line number; starts from 1"""
        self.column = column
        """column number; starts from 1"""

    def __lt__(self, other):
        if self.line < other.line:
            return True
        elif self.line == other.line:
            return self.column < other.column
        return False

    def __str__(self):
        return "[{line}:{column}]".format(line=self.line, column=self.column)

    @staticmethod
    def from_mark(mark):
        """Returns a `Position` from YAML mark."""
        return Position(mark.line + 1, mark.column + 1)

    @staticmethod
    def from_document_end(document):
        """Returns the last `Position` in the document."""
        lines = document.splitlines()
        line = len(lines)
        column = len(lines[line-1])
        return Position(line, column + 1)

    @staticmethod
    def from_yaml_error(yaml_error):
        """Returns a `Position` from `MarkedYAMLError`."""
        if yaml_error.problem_mark is not None:
            return Position.from_mark(yaml_error.problem_mark)
        else:
            return Position.from_mark(yaml_error.context_mark)


class Span(ComparableMixin):
    """Borders a part of text."""

    def __init__(self, start=None, end=None):
        self.start = start
        """:class:`.Position` indicates the start of the section; inclusive"""
        self.end = end
        """:class:`.Position` indicates the end of the section; exclusive"""

    def __str__(self):
        return "{start}-{end}".format(
            start=self.start,
            end=self.end
        )

    def __lt__(self, other):
        if self.start is None:
            return True
        if other is None or other.start is None:
            return False
        return self.start < other.start

    @staticmethod
    def from_event(event):
        """Constructs `Span` from YAML `event`."""
        return Span.from_marks(event.start_mark, event.end_mark)

    @staticmethod
    def from_marks(start_mark, end_mark):
        """Constructs `Span` from YAML marks."""
        start = Position.from_mark(start_mark)
        end = Position.from_mark(end_mark)
        return Span(start, end)
