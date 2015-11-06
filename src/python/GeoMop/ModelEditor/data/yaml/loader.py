"""
Package for generating DataNode structure from YAML document.
"""

__author__ = 'Tomas Krizek'

import re
import copy

import yaml as pyyaml

from .constructor import construct_scalar
from .resolver import resolve_scalar_tag
from ..util import TextValue
from helpers import Position, Span
from ..data_node import ScalarNode, CompositeNode, NodeOrigin
from helpers import Notification, NotificationHandler


class Loader:
    """Generates DataNode structure from YAML document."""
    def __init__(self, notification_handler=None):
        """Initializes the loader with NotificationHandler."""
        self._event = None
        self._event_generator = iter([])
        self._document = None
        self._iterate_events = True
        self._fatal_error_node = None
        self.anchors = {}
        self.notification_handler = notification_handler
        if self.notification_handler is None:
            self.notification_handler = NotificationHandler()

    def load(self, document):
        """Loads the YAML document and returns the root DataNode."""
        if document is None:
            return None
        self._fatal_error_node = None
        self._iterate_events = True
        self.anchors = {}
        self._document = document
        self._event_generator = pyyaml.parse(self._document)

        self._next_parse_event()  # StreamStartEvent
        self._next_parse_event()  # DocumentStartEvent
        self._next_parse_event()  # first actual event

        root = self._create_node()
        if self._fatal_error_node:
            if isinstance(root, CompositeNode):
                # pylint: disable=no-member
                root.set_child(self._fatal_error_node)
                root.span.end = self._fatal_error_node.end
        return root

    def _next_parse_event(self):
        """
        Attempts to get next parsing event and handles errors. When there
        are no more valid parsing events, event is set to None.
        """
        if self._iterate_events:
            try:
                # get next parsing event
                self._event = next(self._event_generator)
            except pyyaml.MarkedYAMLError as error:
                # handle parsing error
                self._event = None
                self._iterate_events = False
                start_pos = Position.from_yaml_error(error)
                self._create_fatal_error_node(start_pos)
                notification = Notification.from_name('SyntaxFatalError')
                notification.span = self._fatal_error_node.span
                self.notification_handler.report(notification)
            except StopIteration:
                # handle end of parsing events
                self._event = None
                self._iterate_events = False

    def _create_node(self, parent=None):
        """Create a DataNode from parsing events."""
        anchor = self._extract_anchor()
        tag = self._extract_tag()
        node = None
        if tag is not None:
            node = self._create_node_by_tag(tag)
        elif isinstance(self._event, pyyaml.MappingStartEvent):
            node = self._create_record_node()
        elif isinstance(self._event, pyyaml.SequenceStartEvent):
            node = self._create_array_node()
        elif isinstance(self._event, pyyaml.ScalarEvent):
            node = self._create_scalar_node()
        elif isinstance(self._event, pyyaml.AliasEvent):
            node = self._create_alias_node(anchor)
        if node is not None:
            node.parent = parent
            # register anchor if it exists and node is not an alias
            if anchor is not None and node.ref is None:
                self._register_anchor(anchor, node)
        return node

    def _extract_anchor(self):
        """Extracts `TextValue` of anchor from the current event."""
        value = getattr(self._event, 'anchor', None)
        if value is None or value in ['*', '&']:
            return None
        anchor = TextValue(value)
        symbol = '&'
        if isinstance(self._event, pyyaml.AliasEvent):
            symbol = '*'
        anchor.span = self._extract_property_span(self._event.start_mark, symbol)
        return anchor

    def _extract_tag(self):
        """Extracts `TextValue` of tag from the current event."""
        value = getattr(self._event, 'tag', None)
        if value is None or value == '!':
            return None
        tag = TextValue(value)
        tag.span = self._extract_property_span(self._event.start_mark, '!')
        return tag

    def _extract_property_span(self, start_mark, symbol):
        """
        Create a `Span` from the first `symbol` at `start_mark`.
        Used to get the span of node properties like anchors or tags.
        """
        # set document to start at start_mark
        lines = self._document.splitlines()
        line_index = start_mark.line
        line = lines[line_index]
        line = line[start_mark.column:]  # first line offset
        expr = '[{symbol}]([a-zA-Z0-9_:-]+)'.format(symbol=symbol)
        regex = re.compile(expr)

        match = regex.search(line)
        if match is not None:  # set correct offset of match
            start_column = start_mark.column + match.start(1)
            end_column = start_mark.column + match.end(1)

        while match is None and line_index < len(lines) - 1:
            line_index += 1
            line = lines[line_index]
            match = regex.search(line)

        start_column = locals().get('start_column', match.start(1))
        end_column = locals().get('end_column', match.end(1))
        start = Position(line_index + 1, start_column + 1)
        end = Position(line_index + 1, end_column + 1)
        return Span(start, end)

    def _create_node_by_tag(self, tag):
        """
        Creates either an abstract record (for app specific tags: !) or
        scalar value of specified type (yaml tags - !!, tag:yaml.org,2002:)
        """
        if tag.value.startswith('tag:yaml.org,2002:'):
            node = self._create_scalar_node()
        else:  # abstract record
            node = self._create_abstract_record(tag)
        return node

    def _create_scalar_node(self):
        """Creates a ScalarNode."""
        node = ScalarNode()
        tag = self._event.tag
        if tag is None or not tag.startswith('tag:yaml.org,2002:'):
            tag = resolve_scalar_tag(self._event.value)
        node.span = Span.from_event(self._event)
        try:
            node.value = construct_scalar(self._event.value, tag)
        except Exception as error:
            notification = Notification.from_name('ConstructScalarError', error.args[0])
            notification.span = node.span
            self.notification_handler.report(notification)
            return node
        if node.value is None:
            # alter position of empty node (so it can be selected)
            node.span.end.column += 1
        return node

    def _create_abstract_record(self, tag):
        """Creates abstract record from parsing events."""
        tag.value = tag.value[1:]  # remove leading !
        if isinstance(self._event, pyyaml.MappingStartEvent):
            # classic abstract record node
            node = self._create_record_node()
        elif isinstance(self._event, pyyaml.ScalarEvent):
            # scalar abstract node - either empty, or for autoconversion
            temp_node = self._create_scalar_node()
            if temp_node.value is None:
                # empty node - construct as mapping
                node = CompositeNode(True)
                node.span = temp_node.span
            else:  # may be used for autoconversion
                node = temp_node
        elif isinstance(self._event, pyyaml.SequenceStartEvent):
            node = self._create_array_node()
        node.type = tag
        return node

    def _create_record_node(self):
        """Creates a record node."""
        node = CompositeNode(True)
        start_mark = self._event.start_mark
        end_mark = self._event.end_mark
        self._next_parse_event()
        # create children
        while (self._event is not None and
               not isinstance(self._event, pyyaml.MappingEndEvent)):
            key = self._create_record_key()
            self._next_parse_event()  # value event
            if not key:  # if key is invalid
                continue
            elif key.value == '<<':  # handle merge
                self._perform_merge(key, node)
                self._next_parse_event()
                continue
            if self._event is None:
                break  # something went wrong, abandon ship!
            child_node = self._create_node(node)
            self._next_parse_event()
            if child_node is None:  # i.e. unresolved alias
                continue
            child_node.key = key
            node.set_child(child_node)
        if self._event is not None:  # update end_mark when map ends correctly
            end_mark = self._event.end_mark
        elif node.children:
            end_mark = node.children[-1].span.end
            end_mark.line -= 1
            end_mark.column -= 1
        node.span = Span.from_marks(start_mark, end_mark)
        return node

    def _create_record_key(self):
        """Creates `TextValue` of record key."""
        # check if key is scalar
        if not isinstance(self._event, pyyaml.ScalarEvent):
            start_pos = Position.from_mark(self._event.start_mark)
            self._create_fatal_error_node(start_pos)
            notification = Notification.from_name('ComplexRecordKey')
            notification.span = self._fatal_error_node.span
            self.notification_handler.report(notification)
            self._event = None
            self._iterate_events = False
            return None
        key = TextValue()
        key.value = self._event.value
        key.span = Span.from_event(self._event)
        return key

    def _perform_merge(self, key, node):
        """Performs merge operation on record node."""
        if isinstance(self._event, pyyaml.SequenceStartEvent):
            self._next_parse_event()
            while (self._event is not None and
                   not isinstance(self._event, pyyaml.SequenceEndEvent)):
                self._merge_into_node(node)
                self._next_parse_event()
        else:
            self._merge_into_node(node)

    def _merge_into_node(self, node):
        """Merges a single alias node into this record node."""
        # allow only alias nodes to be merged
        if not isinstance(self._event, pyyaml.AliasEvent):
            temp = self._create_node(node)  # skip the node to avoid parsing error
            notification = Notification.from_name('NoReferenceToMerge')
            notification.span = temp.span
            self.notification_handler.report(notification)
            return

        anchor = self._extract_anchor()
        if anchor.value not in self.anchors:
            notification = Notification.from_name('UndefinedAnchor', anchor.value)
            notification.span = anchor.span
            self.notification_handler.report(notification)
            return

        anchor_node = self.anchors[anchor.value]
        # check if anchor_node is a record (mapping)
        not_composite = not isinstance(anchor_node, CompositeNode)
        if not_composite or not anchor_node.explicit_keys:
            notification = Notification.from_name('InvalidReferenceForMerge', anchor.value)
            notification.span = anchor.span
            self.notification_handler.report(notification)
            return

        for child in anchor_node.children:
            if child.key.value not in node.children_keys:
                node.set_child(copy.deepcopy(child))
        return

    def _create_array_node(self):
        """Creates an array node."""
        node = CompositeNode(False)
        start_mark = self._event.start_mark
        end_mark = self._event.end_mark
        self._next_parse_event()
        while (self._event is not None and
               not isinstance(self._event, pyyaml.SequenceEndEvent)):
            key = TextValue(str(len(node.children)))
            child_node = self._create_node(node)
            self._next_parse_event()
            if child_node is None:  # i.e. unresolved alias
                continue
            child_node.key = key
            node.children.append(child_node)
        if self._event is not None:  # update end_mark when array ends correctly
            end_mark = self._event.end_mark
        elif node.children:
            end_mark = node.children[-1].span.end
            end_mark.line -= 1
            end_mark.column -= 1
        node.span = Span.from_marks(start_mark, end_mark)
        return node

    def _create_alias_node(self, anchor):
        """Creates an alias node."""
        if anchor.value not in self.anchors:
            notification = Notification.from_name('UndefinedAnchor', anchor.value)
            notification.span = anchor.span
            self.notification_handler.report(notification)
            return None
        ref = self.anchors[anchor.value]
        node = copy.deepcopy(ref)
        node.anchor = anchor
        node.ref = ref

        # set correct node.span
        node.span = copy.deepcopy(node.anchor.span)
        start = node.span.start
        node.span.start = Position(start.line, start.column - 1)
        return node

    def _register_anchor(self, anchor, node):
        """Registers an anchor to a node."""
        node.anchor = anchor
        if anchor.value in self.anchors:
            notification = Notification.from_name('OverridingAnchor', anchor.value)
            notification.span = anchor.span
            self.notification_handler.report(notification)
        self.anchors[anchor.value] = node

    def _create_fatal_error_node(self, start_pos):
        """
        Creates a non-existing node in the data tree
        to wrap the content of the error (span) in a node.
        """
        node = ScalarNode()
        end_pos = Position.from_document_end(self._document)
        node.span = Span(start_pos, end_pos)
        node.key = TextValue('fatal_error')
        node.origin = NodeOrigin.error
        node.hidden = True
        self._fatal_error_node = node
