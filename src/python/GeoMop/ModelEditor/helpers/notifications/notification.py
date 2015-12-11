"""
Notification related classes.
"""

__author__ = 'Tomas Krizek'

from enum import Enum

from util.locators import Span, Position
from .list import NOTIFICATIONS_BY_CODE, NOTIFICATIONS_BY_NAME


class Notification(Exception):
    """Notification message."""

    @staticmethod
    def from_code(code, *msg_args):
        """Constructs a `Notification` from code and arguments."""
        name = NOTIFICATIONS_BY_CODE[code]['name']
        message = NOTIFICATIONS_BY_CODE[code]['message'].format(*msg_args)
        return Notification(code, name, message)

    @staticmethod
    def from_name(name, *msg_args):
        """Constructs a `Notification` from name and arguments."""
        code = NOTIFICATIONS_BY_NAME[name]['code']
        message = NOTIFICATIONS_BY_NAME[name]['message'].format(*msg_args)
        return Notification(code, name, message)

    class Severity(Enum):
        """Severity of a notification."""
        info = 0
        warning = 1
        error = 2
        fatal = 3

    def __init__(self, code, name, message, span=None):
        """Initializes the notification."""
        super(Notification, self).__init__()
        self.code = code
        self.name = name
        self.message = message
        self._span = Span(Position(1, 1), Position(1, 1))
        self.span = span

    @property
    def title(self):
        """Title of the `Notification`."""
        return self.name

    @property
    def span(self):
        """Span of the `Notification`."""
        return self._span

    @span.setter
    def span(self, span):
        """Sets the span property."""
        if isinstance(span, Span):
            self._span = span

    @property
    def severity(self):
        """Severity of the `Notification`."""
        if 100 <= self.code < 300:
            return Notification.Severity.fatal
        elif self.code < 600:
            return Notification.Severity.error
        elif self.code < 800:
            return Notification.Severity.warning
        else:
            return Notification.Severity.info

    def __str__(self):
        text = "{code:3d} ({name}): {message}"
        return text.format(
            code=self.code,
            name=self.name,
            message=self.message
        )
