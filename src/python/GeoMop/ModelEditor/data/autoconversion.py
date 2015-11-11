"""
GeoMop model auto-conversion module

Ensures auto-conversion of data for specified format.
"""

__author__ = 'Tomas Krizek'

from copy import deepcopy

from .data_node import CompositeNode, NodeOrigin
from .util import TextValue


def autoconvert(node, input_type):
    """
    Performs recursive auto-conversion on root node.

    Auto-correction:
        1. If Array is expected and scalar/record is found, encapsulate it
           in Array(s).
        2. If Record is expected and scalar/array is found, check
           reducible_to_key. If present, create the Record.
        3. If AbstractRecord is expected and scalar/array is found, check if
           default_descendant fits rule 2.
    """
    root = deepcopy(node)
    _autoconvert_crawl(root, input_type)
    return root


def _autoconvert_crawl(node, input_type):
    """
    Recursively crawls through the tree structure and tries to auto-convert
    values to the expected type.
    """
    if input_type is None:
        return

    if input_type['base_type'] == 'AbstractRecord':
        try:
            it_concrete = input_type['implementations'][node.type.value]
        except (KeyError, AttributeError):
            try:
                it_concrete = input_type['default_descendant']
            except KeyError:
                return
        _autoconvert_crawl(node, it_concrete)
    elif input_type['base_type'] == 'Array':
        if not isinstance(node, CompositeNode):
            return
        children = list(node.children)
        node.children.clear()
        for child in children:
            ac_child = _get_autoconverted(child, input_type['subtype'])
            node.set_child(ac_child)
            _autoconvert_crawl(ac_child, input_type['subtype'])
    elif input_type['base_type'] == 'Record':
        if not isinstance(node, CompositeNode):
            return
        children = list(node.children)
        node.children.clear()
        for child in children:
            try:
                child_it = input_type['keys'][child.key.value]['type']
            except (KeyError, AttributeError):
                node.set_child(child)
                continue
            else:
                ac_child = _get_autoconverted(child, child_it)
                node.set_child(ac_child)
                _autoconvert_crawl(ac_child, child_it)
    return


def _get_autoconverted(node, input_type):
    """
    Auto-conversion of array and record types.

    Arrays are expanded to the expected dimension.
    Records are initialized from the reducible_to_key value.
    """
    if input_type is None:
        return node
    is_array = isinstance(node, CompositeNode) and not node.explicit_keys
    is_record = isinstance(node, CompositeNode) and node.explicit_keys
    if input_type['base_type'] == 'Array' and not is_array:
        dim = _get_expected_array_dimension(input_type)
        return _expand_value_to_array(node, dim)
    elif input_type['base_type'].endswith('Record') and not is_record:
        return _expand_reducible_to_key(node, input_type)
    else:
        return node


def _get_expected_array_dimension(input_type):
    """Returns the expected dimension of the input array."""
    dim = 0
    while input_type['base_type'] == 'Array':
        dim += 1
        input_type = input_type['subtype']
    return dim


def _expand_value_to_array(node, dim):
    """Expands node value to specified dimension."""
    while dim > 0:
        array_node = CompositeNode(False, node.key, node.parent)
        array_node.span = node.span
        node.parent = array_node
        node.key = TextValue('0')
        if node.input_type is not None:
            array_node.input_type = node.input_type
            node.input_type = array_node.input_type['subtype']
        array_node.children.append(node)
        array_node.origin = NodeOrigin.ac_array
        node = array_node
        dim -= 1
    return node


def _expand_reducible_to_key(node, input_type):
    """Initializes a record from the reducible_to_key value."""
    if input_type is None:
        return
    try:
        key = input_type['default_descendant']['reducible_to_key']
        child_input_type = input_type['default_descendant']
    except (KeyError, TypeError):
        try:
            key = input_type['reducible_to_key']
            child_input_type = input_type
        except KeyError:
            return node

    if key is None:
        return node

    record_node = CompositeNode(True, node.key, node.parent)
    record_node.span = node.span
    if hasattr(node, 'type'):
        record_node.type = node.type
        node.type = None
    node.parent = record_node
    node.origin = NodeOrigin.ac_reducible_to_key
    node.key = TextValue(key)
    if node.input_type is not None:
        record_node.input_type = node.input_type
        node.input_type = child_input_type['keys'][key]['type']
    record_node.children.append(node)
    return record_node
