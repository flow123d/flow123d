"""
Data Node package

Contains classes for representing the tree structure of config files.
"""

__author__ = 'Tomas Krizek'

from enum import Enum

from .util import TextValue

DEBUG_MODE = True
"""changes the behaviour to debug mode"""


class DataNode:
    """
    Represents a node in the tree structure.

    The complete tree is represented by its root node.
    """
    def __init__(self, key=None, parent=None):
        self.ref = None
        """reference to another node"""
        self.parent = parent
        """parent node"""
        self.key = key
        """key (name) of this node"""
        if self.key is None:
            self.key = TextValue()
        self.input_type = None
        """input type specified by format"""
        self.span = None
        """borders the position of this node in input text"""
        self.anchor = None
        """anchor of the node `TextValue`"""
        self.is_flow = False
        """outer structure of the node is in flow format"""
        self.delimiters = None
        """Outer border of node, span start point to first delimiter and to second"""
        self.origin = NodeOrigin.structure
        """indicates how node was created"""
        self.hidden = False
        """whether node is hidden in tree structure"""
        self._options = []

    @property
    def absolute_path(self):
        """the absolute path to this node"""
        return self.generate_absolute_path()

    def generate_absolute_path(self, descendant_path=None):
        """generates absolute path to this node, used recursively"""
        if self.parent is None:
            if descendant_path is None:
                descendant_path = ""
            return "/" + descendant_path
        if descendant_path is None:
            path = str(self.key.value)
        else:
            path = str(self.key.value) + "/" + descendant_path
        return self.parent.generate_absolute_path(path)

    @property
    def options(self):
        """possible options to hint in autocomplete"""
        options = self._options
        # if DEBUG_MODE and self.span is not None:
        #     # return start, end position as options
        #     options = ["start: {start}".format(start=self.span.start),
        #                "end: {end}".format(end=self.span.end)]
        return options

    @options.setter
    def options(self, options):
        """Autocomplete options setter."""
        self._options = options

    def get_node_at_position(self, position):
        """Retrieves DataNode at specified position."""
        raise NotImplementedError

    def is_child_on_line(self, line):
        """
        Return if in set line is some child
        """
        return False

    def is_jsonchild_on_line(self, line):
        """
        Return if in set line is some json child
        """
        return self.is_flow

    def get_node_at_path(self, path):
        """returns node at given path"""
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
            if not isinstance(node, CompositeNode) or node.get_child(key) is None:
                raise LookupError("Node {key} does not exist in {location}"
                                  .format(key=key, location=node.absolute_path))
            node = node.get_child(key)
        return node

    def get_info_text_data(self, is_parent=False):
        """
        Returns data necessary to generate info_text.

        `is_parent` should be set to True when the generated info_text should be for this node,
        instead of its parent node
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
            while (node.origin == NodeOrigin.ac_array or
                   not (isinstance(node, CompositeNode) and node.explicit_keys is True) or
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

    @property
    def start(self):
        """start of node, including its key"""
        start = self.span.start
        if self.key is not None and self.key.span is not None:
            start = self.key.span.start
        return start

    @property
    def end(self):
        """Returns the end of this node."""
        return self.span.end


class CompositeNode(DataNode):
    """Represents a composite node in the tree structure."""
    def __init__(self, explicit_keys, key=None, parent=None):
        super(CompositeNode, self).__init__(key, parent)
        self.children = []
        """list of children nodes"""
        self.explicit_keys = explicit_keys
        """boolean; indicates whether keys are specified (record) or
        implicit (array)"""
        self.type = None
        """specifies the type of AbstractRecord"""

    def get_node_at_position(self, position):
        """Retrieves DataNode at specified position."""
        node = None
        if self.start <= position <= self.end:
            node = self
            for child in self.children:
                descendant = child.get_node_at_position(position)
                if descendant is not None:
                    node = descendant
                    break
        return node

    def __str__(self):
        text = super(CompositeNode, self).__str__()
        children_keys = [str(child.key.value) for child in self.children]
        text += "  children_keys: {children_keys}\n".format(
            children_keys=', '.join(children_keys)
        )
        return text

    def get_child(self, key):
        """
        Returns a child node for the given key; None if key doesn't
        exist.
        """
        for child in self.children:
            if key == child.key.value:
                return child
        return None

    def set_child(self, node):
        """
        Sets the specified node as child of this node. If the key already
        exists, the other child node is replaced by this child_node.
        """
        node.parent = self

        for i, child in enumerate(self.children):
            if child.key.value == node.key.value:
                self.children[i] = node
                return
        # still not ended - new key
        self.children.append(node)

    def is_child_on_line(self, line):
        """
        Return if in set line is some child
        """
        for child in self.children:
            if child.start.line <= line and child.end.line >= line:
                return True
        return False

    def is_jsonchild_on_line(self, line):
        """
        Return if in set line is some json child
        """
        for child in self.children:
            if child.start.line <= line and child.end.line >= line:
                if child.is_flow:
                    return True
                return child.is_jsonchild_on_line(line)
        return False

    @property
    def children_keys(self):
        """Returns all children keys."""
        return [child.key.value for child in self.children]

    @property
    def visible_children(self):
        """Returns a list of all visible children nodes."""
        return [child for child in self.children if child.hidden is False]


class ScalarNode(DataNode):
    """Represents a scalar node in the tree structure."""
    def __init__(self, key=None, parent=None, value=None):
        super(ScalarNode, self).__init__(key, parent)
        self.value = value
        """the scalar value"""

    def get_node_at_position(self, position):
        """Retrieves DataNode at specified position."""
        if self.start <= position <= self.end:
            return self
        return None

    def __str__(self):
        text = super(ScalarNode, self).__str__()
        text += "  value: {value}".format(
            value=self.value
        )
        return text


class NodeOrigin(Enum):
    """The origin of data node."""
    structure = 1
    ac_array = 2
    ac_reducible_to_key = 3
    error = 4
