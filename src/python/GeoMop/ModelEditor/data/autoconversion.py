"""
GeoMop model auto-conversion module

Ensures auto-conversion of data for specified format.

.. codeauthor:: Tomas Krizek <tomas.krizek1@tul.cz>
"""

from copy import deepcopy

from helpers import notification_handler, Notification
from util.util import TextValue

from .data_node import DataNode, MappingDataNode, SequenceDataNode
from .format import SCALAR


class AutoConverter:
    """Handle autoconverting layer of data."""

    @staticmethod
    def autoconvert(node, input_type):
        """
        Performs recursive auto-conversion on root node.

        Auto-conversions:
            1. If Array is expected and scalar/record is found, either perform
               a transposition or encapsulate it in Array(s).
            2. If Record is expected and scalar/array is found, check
               reducible_to_key. If present, create the Record.
            3. If AbstractRecord is expected and scalar/array is found, check if
               default_descendant fits rule 2.

        The function also converts data to the expected data types (if possible).
        """
        root = deepcopy(node)
        AutoConverter._autoconvert_crawl(root, input_type)
        return root

    @staticmethod
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
            AutoConverter._autoconvert_crawl(node, it_concrete)
        elif input_type['base_type'] == 'Array':
            if node.implementation != DataNode.Implementation.sequence:
                return
            children = list(node.children)
            node.children.clear()
            for child in children:
                ac_child = AutoConverter._get_autoconverted(child, input_type['subtype'])
                node.set_child(ac_child)
                AutoConverter._autoconvert_crawl(ac_child, input_type['subtype'])
        elif input_type['base_type'] == 'Record':
            if node.implementation != DataNode.Implementation.mapping:
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
                    ac_child = AutoConverter._get_autoconverted(child, child_it)
                    node.set_child(ac_child)
                    AutoConverter._autoconvert_crawl(ac_child, child_it)
        elif input_type['base_type'] in SCALAR:
            ScalarConverter.convert(node, input_type)

        return

    @staticmethod
    def _get_autoconverted(node, input_type):
        """
        Auto-conversion of array and record types.

        Arrays are expanded to the expected dimension.
        Records are initialized from the reducible_to_key value.
        """
        if input_type is None:
            return node
        is_array = node.implementation == DataNode.Implementation.sequence
        is_record = node.implementation == DataNode.Implementation.mapping
        if input_type['base_type'] == 'Array' and not is_array:
            return transposer.make_transposition(node, input_type)
        elif input_type['base_type'].endswith('Record') and not is_record:
            return AutoConverter._expand_reducible_to_key(node, input_type)
        else:
            return node

    @staticmethod
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

        record_node = MappingDataNode(node.key, node.parent)
        record_node.span = node.span
        if hasattr(node, 'type'):
            record_node.type = node.type
            node.type = None
        node.parent = record_node
        node.origin = DataNode.Origin.ac_reducible_to_key
        node.key = TextValue(key)
        if node.input_type is not None:
            record_node.input_type = node.input_type
            node.input_type = child_input_type['keys'][key]['type']
        record_node.children.append(node)
        return record_node


class ScalarConverter:
    """Convert scalar values to their expected types."""

    @staticmethod
    def convert(node, input_type):
        """Convert value of Scalar node to expected type.

        node: :py:class:`DataNode` data structure
        input_type: definition of input_type
        """
        conversions = {
            'Bool': ScalarConverter._convert_to_bool,
            'Integer': ScalarConverter._convert_to_int,
            'Double': ScalarConverter._convert_to_float,
            'String': ScalarConverter._convert_to_string,
            'FileName': ScalarConverter._convert_to_string,
            'Selection': ScalarConverter._convert_to_string,
        }

        base_type = input_type['base_type']
        if base_type in conversions and node.value is not None:
            try:
                value = conversions[base_type](node.value)
            except ValueError:
                notification = Notification.from_name('ValueConversionError', node.value, base_type)
                notification.span = node.span
                notification_handler.report(notification)
                return
            node.value = value

    @staticmethod
    def _convert_to_bool(value):
        """Convert given value to bool."""
        if not isinstance(value, bool):
            return value.lower() in ("true", "1")
        return value

    @staticmethod
    def _convert_to_int(value):
        """Convert given value to int.

        :raises: ValueError - if the value can not be converted to integer
        """
        if not isinstance(value, int):
            return int(ScalarConverter._convert_to_float(value))
        return value

    @staticmethod
    def _convert_to_float(value):
        """Convert given value to float.

        :raises: ValueError - if the value can not be converted to float
        """
        if not isinstance(value, float):
            return float(value)
        return value

    @staticmethod
    def _convert_to_string(value):
        """Convert the given value to string."""
        if not isinstance(value, str):
            return str(value)
        return value


class Transposer:
    """Handle the transposition autoconversion.

    This conversion happens when Array is expected, but Scalar od Record is encountered.
    Scalar values are simply encapsulated in an array of the correct dimension.

    If there is a Record, check if the sizes of all the unexpected arrays inside it
    (recursively) match. If so, proceed with the conversion. Transpose the Record into an
    Array of Records. The size of an Array is determined by the previous step. Finally,
    ensure all of the unexpected Arrays are replaced by a single value at the correct
    position.
    """

    paths_to_convert = None
    """a list of path to keys which are to be converted"""
    array_size = None
    """the size of the array tobe created by transposition"""
    current_path = None
    """list of keys that lead to the current node"""

    @classmethod
    def init(cls):
        """Initialize class for operation."""
        cls.paths_to_convert = []
        cls.array_size = None
        cls.current_path = ['.']

    @classmethod
    def make_transposition(cls, node, input_type):
        """Transpose a record or scalar into an array."""
        assert input_type['base_type'] == 'Array', "Only Array can be a result of transposition"
        cls.init()

        # if node is scalar, convert it to array
        if node.implementation == DataNode.Implementation.scalar:
            return cls._expand_value_to_array(node)

        # verify that subtype is record
        subtype = input_type['subtype']
        if subtype['base_type'] != 'Record':
            notification = Notification.from_name('UnsupportedTransposition',
                                                  input_type['base_type'])
            notification.span = node.span
            notification_handler.report(notification)
            return node
        assert node.implementation == DataNode.Implementation.mapping,\
            "Can not perform transposition on array"

        # get array size
        try:
            cls._get_transformation_array_size(node, subtype)
        except Notification as notification:
            notification_handler.report(notification)
            return node
        if cls.array_size is None:
            cls.array_size = 1

        # create array
        array_node = SequenceDataNode(node.key, node.parent)
        array_node.span = node.span
        array_node.input_type = node.input_type
        array_node.origin = DataNode.Origin.ac_array
        template_node = deepcopy(node)
        template_node.parent = array_node
        template_node.input_type = subtype
        template_node.origin = DataNode.Origin.ac_transposition

        # create and transpose items of the array
        for i in range(cls.array_size):
            child_node = deepcopy(template_node)
            child_node.key = TextValue(str(i))
            # convert array to value
            for path in cls.paths_to_convert:
                node_to_convert = child_node.get_node_at_path(path)
                converted_node = node_to_convert.children[i]
                converted_node.parent = node_to_convert.parent
                converted_node.key = node_to_convert.key
                node_to_convert.parent.set_child(converted_node)
            array_node.children.append(child_node)

        return array_node

    @classmethod
    def _get_transformation_array_size(cls, node, input_type):
        """Return transformation array size."""
        # find a children node that has an array instead of record or scalar
        for child in node.children:
            # the key is not specified in input type
            if 'keys' not in input_type or child.key.value not in input_type['keys']:
                continue

            child_type = input_type['keys'][child.key.value]['type']

            if child.implementation == DataNode.Implementation.sequence:
                if child_type['base_type'] == 'Record':
                    notification = Notification.from_name("InvalidTransposition")
                    notification.span = child.span
                    raise notification
                elif child_type['base_type'] != 'Array':
                    if cls.array_size is None:
                        cls.array_size = len(child.children)
                        cls.paths_to_convert.append('/'.join(cls.current_path + [child.key.value]))
                    elif cls.array_size != len(child.children):
                        notification = Notification.from_name(
                            "DifferentArrayLengthForTransposition")
                        notification.span = child.span
                        raise notification
                    else:
                        cls.paths_to_convert.append('/'.join(cls.current_path + [child.key.value]))

            # verify array size recursively
            cls.current_path.append(child.key.value)
            cls._get_transformation_array_size(child, child_type)
            cls.current_path.pop()

        return cls.array_size

    @staticmethod
    def _expand_value_to_array(node):
        """Expands node value to an array."""
        array_node = SequenceDataNode(node.key, node.parent)
        array_node.span = node.span
        node.parent = array_node
        node.key = TextValue('0')
        if node.input_type is not None:
            array_node.input_type = node.input_type
            node.input_type = array_node.input_type['subtype']
        array_node.children.append(node)
        array_node.origin = DataNode.Origin.ac_array
        return array_node


# initialize module
autoconvert = AutoConverter.autoconvert
transposer = Transposer()

__all__ = ['autoconvert', 'AutoConverter']
