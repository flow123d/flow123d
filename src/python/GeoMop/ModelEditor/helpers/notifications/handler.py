"""
Notification Handler module for capturing and reporting notifications
that occur when processing the data structure.

.. codeauthor:: Tomas Krizek <tomas.krizek1@tul.cz>
"""


class NotificationHandler:
    """
    Handles notifications that can occur during parsing, loading or validation
    of data structure.
    """

    def __init__(self):
        self._notifications = []

    def clear(self):
        """clears the notification buffer"""
        self._notifications = []

    @property
    def notifications(self):
        """Sorted array of notifications."""
        return sorted(self._notifications, key=lambda notification: (
            notification.span, -notification.severity.value))

    def report(self, notification):
        """Reports a notification."""
        if notification.title == 'SilencedNotification':
            return
        self._notifications.append(notification)


notification_handler = NotificationHandler()
