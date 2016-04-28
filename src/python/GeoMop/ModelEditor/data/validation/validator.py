"""Validator for Flow123D data structure

.. codeauthor:: Tomas Krizek <tomas.krizek1@tul.cz>
"""
from helpers import Notification
from util import TextValue, Span

from geomop_project import Parameter

from . import checks
from ..data_node import DataNode
from ..format import is_scalar, is_param


class Validator:
    """Handles data structure validation."""

    def __init__(self, notification_handler):
        """Initializes the validator with a NotificationHandler."""
        self.notification_handler = notification_handler
        self.valid = True
        self.params = []

    def validate(self, node, input_type):
        """
        Performs data validation of node with the specified input_type.

        Validation is performed recursively on all children nodes as well.
        Options are added to nodes where applicable (record keys, selection, ...).

        Returns True when all data was correctly validated, False otherwise.
        Attribute errors contains a list of occurred errors.
        """
        self.valid = True
        self.params = []
        self._validate_node(node, input_type)
        return self.valid

    def _validate_node(self, node, input_type):
        """
        Determines if node contains correct value.

        Method verifies node recursively. All descendant nodes are checked.
        """
        if node is None:
            raise Notification.from_name('ValidationError', 'Invalid node (None)')

        # parameters
        # TODO: enable parameters in unknown IST?
        if hasattr(node, 'value'):
            match = is_param(node.value)
            if match:
                # extract parameters
                new_param = Parameter(match.group(1))
                exists = False
                for param in self.params:
                    if param.name == new_param.name:
                        exists = True
                        break
                if not exists:
                    self.params.append(new_param)
                node.input_type = input_type
                # assume parameters are correct, do not validate further
                return

        if input_type['base_type'] != 'Abstract' and hasattr(node, 'type') \
                and node.type is not None and 'implemented_abstract_record' not in input_type:
            notification = Notification.from_name('UselessTag', node.type.value)
            notification.span = node.type.span
            self.notification_handler.report(notification)

        node.input_type = input_type
        if is_scalar(input_type):
            self._validate_scalar(node, input_type)
        elif input_type['base_type'] == 'Record':
            self._validate_record(node, input_type)
        elif input_type['base_type'] == 'Abstract':
            self._validate_abstract(node, input_type)
        elif input_type['base_type'] == 'Array':
            self._validate_array(node, input_type)
        else:
            notification = Notification.from_name('InputTypeNotSupported',
                                                  input_type['base_type'])
            self._report_notification(notification)

    def _validate_scalar(self, node, input_type):
        """Validates a Scalar node."""
        if input_type['base_type'] == 'Selection':
            node.options = input_type['values']
        try:
            checks.check_scalar(node, input_type)
        except Notification as notification:
            if notification.name in ['InvalidSelectionOption', 'ValueTooBig', 'ValueTooSmall',
                                     'ValidationTypeError']:
                notification.span = node.span
            else:
                notification.span = get_node_key(node).notification_span
            self._report_notification(notification)

    def _validate_record(self, node, input_type):
        """Validates a Record node."""
        if not node.implementation == DataNode.Implementation.mapping:
            notification = Notification.from_name('ValidationTypeError', 'Record')
            notification.span = get_node_key(node).notification_span
            self._report_notification(notification)
            return
        keys = node.children_keys
        node.options = input_type['keys'].keys()
        keys.extend(input_type['keys'].keys())
        for key in set(keys):
            if node.origin == DataNode.Origin.error:
                continue
            child = node.get_child(key)
            try:
                checks.check_record_key(node.children_keys, key, input_type)
            except Notification as notification:
                if notification.name == 'UnknownRecordKey':
                    notification.span = child.notification_span
                else:
                    notification.span = get_node_key(node).notification_span
                self._report_notification(notification)
            else:
                if child is not None:
                    child_input_type = input_type['keys'][key]['type']
                    self._validate_node(child, child_input_type)

    def _validate_abstract(self, node, input_type):
        """Validates an AbtractRecord node."""
        try:
            concrete_type = checks.get_abstractrecord_type(node, input_type)
        except Notification as notification:
            if notification.name == 'InvalidAbstractType':
                notification.span = node.type.span
            else:
                notification.span = get_node_key(node).notification_span
            self._report_notification(notification)
        else:
            if node.type is None:
                # if default_descendant defines the Abstract type, add it to data structure
                node.type = TextValue()
                node.type.value = concrete_type.get('name')
                node.type.span = Span(node.span.start, node.span.start)
            concrete_type['implemented_abstract_record'] = input_type
            node.input_type = concrete_type
            self._validate_record(node, concrete_type)

    def _validate_array(self, node, input_type):
        """Validates an Array node."""
        if not node.implementation == DataNode.Implementation.sequence:
            notification = Notification.from_name('ValidationTypeError', 'Array')
            notification.span = get_node_key(node).notification_span
            self._report_notification(notification)
            return
        try:
            checks.check_array(node.children, input_type)
        except Notification as notification:
            notification.span = get_node_key(node).notification_span
            self._report_notification(notification)
        for child in node.children:
            self._validate_node(child, input_type['subtype'])

    def _report_notification(self, notification):
        """Reports a notification."""
        if notification.severity.value >= Notification.Severity.error.value:
            self.valid = False
        self.notification_handler.report(notification)


def get_node_key(node):
    """Return node that has originated from the text structure (not autoconversion)."""
    while node.origin != DataNode.Origin.structure:
        node = node.parent
    return node
