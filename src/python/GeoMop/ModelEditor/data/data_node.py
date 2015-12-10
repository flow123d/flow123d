"""
Data Node package

Contains classes for representing the tree structure of config files.

.. codeauthor:: Tomas Krizek <tomas.krizek1@tul.cz>
"""

from enum import Enum

from util.util import TextValue


class DataNode:
    """
    Represents a node in the tree structure.

    The complete tree is represented by its root node.
    """

    class Implementation(Enum):
        """Implementation type of :py:class:`DataNode`.

        The different implementations are defined by the three node types in YAML.
        """
        scalar = 1
        sequence = 2
        mapping = 3

    class Origin(Enum):
        """The origin of data node."""
        structure = 1
        ac_array = 2
        ac_reducible_to_key = 3
        error = 4

    class StructureType(Enum):
        """The type of node in the text structure."""
        scalar = 1
        array = 2
        dict = 3
        json_array = 4
        json_dict = 5

    def __init__(self, key=None, parent=None):
        self.ref = None
        """reference to another :py:class:`DataNode`"""
        self.implementation = None
        """the type node implementation - see :py:class:`DataNode.Implementation`"""
        self.parent = parent
        """parent :py:class:`DataNode`"""
        self.children = []
        """list of children nodes"""
        self.key = key
        """key (name) of this node (:py:class:`TextValue <util.util.TextValue>`)"""
        if self.key is None:
            self.key = TextValue()
        self.input_type = None
        """input type specified by format"""
        self.span = None
        """the :py:class:`Span <util.locators.Span>` of this node's value in text"""
        self.anchor = None
        """anchor of the node (:py:class:`TextValue <util.util.TextValue>`)"""
        self.is_flow = False
        """outer structure of the node is in flow format"""
        self.delimiters = None
        """outer border of node, span start point to first delimiter and to second; used for
        flow style nodes"""
        self.origin = DataNode.Origin.structure
        """indicates whether node is written in text or was created by an autoconversion
        (:py:class:`DataNode.Origin`)"""
        self.hidden = False
        """whether node is hidden in tree structure"""
        self.value = None
        """the scalar value of the node (:py:class:`TextValue <util.util.TextValue>`"""
        self.type = None
        """specifies the type of AbstractRecord"""  # TODO is type TextValue?

    @property
    def absolute_path(self):
        """the absolute path to this node"""
        return self.generate_absolute_path()

    @property
    def start(self):
        """the beginning :py:class:<Position <util.locators.Position>` of node
        (including its key)"""
        start = self.span.start
        if self.key is not None and self.key.span is not None:
            start = self.key.span.start
        return start

    @property
    def end(self):
        """the end :py:class:<Position <util.locators.Position>` of node"""
        return self.span.end

    @property
    def children_keys(self):
        """all children node keys"""
        return [child.key.value for child in self.children]

    @property
    def visible_children(self):
        """all visible children nodes"""
        return [child for child in self.children if child.hidden is False]

    @property
    def notification_span(self):
        """span for notification, preferably returns key span, falls back to span"""
        if self.key is not None and self.key.span is not None:
            return self.key.span
        else:
            return self.span

    def generate_absolute_path(self, descendant_path=None):
        """generate absolute path to this node, used recursively"""
        if self.parent is None:
            if descendant_path is None:
                descendant_path = ""
            return "/" + descendant_path
        if descendant_path is None:
            path = str(self.key.value)
        else:
            path = str(self.key.value) + "/" + descendant_path
        return self.parent.generate_absolute_path(path)

    def get_node_at_position(self, position):
        """Retrieve DataNode at specified position."""
        raise NotImplementedError

    def get_child(self, key):
        """Return a child node for the given key.

        :return: child node of given key
        :rtype: :py:class:`DataNode` or ``None``
        """
        for child in self.children:
            if key == child.key.value:
                return child
        return None

    def set_child(self, node):
        """Set the given node as child of this node.

        If the key already exists, replace the original child node.

        :param: :py:class:`DataNode` to be set as a child
        """
        if self.implementation == DataNode.Implementation.scalar:
            raise TypeError("Scalar node can not have children.")

        node.parent = self

        # does the key already exists? if so, replace it
        for i, child in enumerate(self.children):
            if child.key.value == node.key.value:
                self.children[i] = node
                return

        # the key does not exists, create a new child node
        self.children.append(node)
        return

    def is_child_on_line(self, line):
        """Return if in set line is some child"""
        return False

    def is_jsonchild_on_line(self, line):
        """Return if in set line is some json child"""
        return self.is_flow

    def get_node_at_path(self, path):
        """Find node at the given path.

        :param str path: absolute or relative path to the node, i.e. ``/problem/pause_on_run``
        :return: the found node
        :rtype: DataNode
        :raises: LookupError - if no node is found
        """
        # TODO it would make more sense to return None instead of raising LookupError
        # pylint: disable=no-member
        if path is None:
            raise LookupError("No path provided")
        node = self
        if path.startswith(self.absolute_path):  # absolute path
            path = path[len(self.absolute_path):]
        elif path.startswith('/'):  # absolute path with different location
            while node.parent is not None:
                node = node.parent  # crawl up to root

        for key in path.split('/'):
            if not key or key == '.':
                continue
            elif key == '..':
                node = node.parent
                continue
            if node.get_child(key) is None:
                raise LookupError("Node {key} does not exist in {location}"
                                  .format(key=key, location=node.absolute_path))
            node = node.get_child(key)
        return node

    def get_info_text_data(self, is_parent=False):
        """Return data necessary to generate info_text.

        :param bool is_parent: if set to True, generate info_text for this node,
           instead of its parent node
        :return: necessary input_type ids and values to generate info_text
        :rtype: dict
        """
        # pylint: disable=no-member
        abstract_id = None
        selected_item = None
        record_id = None
        selected_key = None

        input_type = 'None' if self.input_type is None else self.input_type.get('base_type')
        if input_type == 'Selection':
            selected_item = self.value

        if is_parent and input_type in ['AbstractRecord', 'Record']:
            node = self
        else:
            # find first parent record node
            prev_node = self
            node = self.parent if self.parent is not None else self
            while (node.origin == DataNode.Origin.ac_array or
                   node.implementation != DataNode.Implementation.mapping or
                   node.input_type is None):
                if node.parent is None:
                    break
                prev_node = node
                node = node.parent
            selected_key = prev_node.key.value

        if node.input_type is not None:
            if 'implemented_abstract_record' in node.input_type:
                abstract_id = node.input_type['implemented_abstract_record']['id']
            if node.input_type['base_type'] == 'AbstractRecord':
                abstract_id = node.input_type['id']
            else:
                record_id = node.input_type['id']

        return {
            'record_id': record_id,
            'selected_key': selected_key,
            'abstract_id': abstract_id,
            'selected_item': selected_item
        }

    def __str__(self):
        text = (
            "{type_} at 0x{address:x}\n"
            "  key: {key}\n"
            "  parent: {parent_type} at 0x{parent_address:x}\n"
            "  ref: {ref}\n"
            "  span: {span}\n"
            "  input_type: {input_type}\n"
        )
        return text.format(
            type_=type(self).__name__,
            address=id(self),
            key=self.key.value,
            parent_type=type(self.parent).__name__,
            parent_address=id(self.parent),
            ref=self.ref,
            span=self.span,
            input_type=self.input_type
        )


class CompositeDataNode(DataNode):
    """Class defines the common behaviour of both Sequence an Mapping nodes."""

    def __str__(self):
        text = super(CompositeDataNode, self).__str__()
        children_keys = [str(child.key.value) for child in self.children]
        text += "  children_keys: {children_keys}\n".format(
            children_keys=', '.join(children_keys)
        )
        return text

    def is_child_on_line(self, line):
        """Return if in set line is some child"""
        for child in self.children:
            if child.start.line <= line <= child.end.line:
                return True
        return False

    def is_jsonchild_on_line(self, line):
        """Return if in set line is some json child"""
        for child in self.children:
            if child.start.line <= line <= child.end.line:
                if child.is_flow:
                    return True
                return child.is_jsonchild_on_line(line)
        return False

    def get_node_at_position(self, position):
        """Retrieve :py:class:`DataNode` at specified
        :py:class:`Position <util.locators.Position>`."""
        node = None
        if self.start <= position <= self.end:
            node = self
            for child in self.children:
                descendant = child.get_node_at_position(position)
                if descendant is not None:
                    node = descendant
                    break
        return node


class MappingDataNode(CompositeDataNode):
    """Implementation of :py:class:`DataNode` for mapping nodes."""

    def __init__(self, key=None, parent=None):
        super(MappingDataNode, self).__init__(key, parent)
        self.implementation = DataNode.Implementation.mapping


class SequenceDataNode(CompositeDataNode):
    """Implementation of :py:class:`DataNode` for sequence nodes."""

    def __init__(self, key=None, parent=None):
        super(SequenceDataNode, self).__init__(key, parent)
        self.implementation = DataNode.Implementation.sequence


class ScalarDataNode(DataNode):
    """Implementation of :py:class:`DataNode` for scalar nodes."""

    def __init__(self, key=None, parent=None, value=None):
        super(ScalarDataNode, self).__init__(key, parent)
        self.implementation = DataNode.Implementation.scalar
        self.value = value

    def get_node_at_position(self, position):
        """Retrieve :py:class:`DataNode` at specified
        :py:class:`Position <util.locators.Position>`."""
        if self.start <= position <= self.end:
            return self
        return None

    def __str__(self):
        text = super(ScalarDataNode, self).__str__()
        text += "  value: {value}".format(
            value=self.value
        )
        return text
