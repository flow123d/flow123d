"""
Basic rules for data validation
"""

__author__ = 'Tomas Krizek'

# pylint: disable=unused-argument

from ..data_node import ScalarNode
import helpers.notifications.notification as ntf


def check_scalar(node, input_type):
    """Checks scalar node value."""
    if not isinstance(node, ScalarNode):
        raise ntf.Notification.from_name('ValidationTypeError', 'Scalar')
    checks = {
        'Integer': check_integer,
        'Double': check_double,
        'Bool': check_bool,
        'String': check_string,
        'Selection': check_selection,
        'FileName': check_filename
    }
    check = checks.get(input_type['base_type'], None)
    if check is None:
        raise ntf.Notification.from_name('ValidationTypeError', 'Scalar')
    return check(node.value, input_type)


def check_integer(value, input_type):
    """Checks if value is an integer within given range."""
    if not isinstance(value, int):
        raise ntf.Notification.from_name('ValidationTypeError', 'Integer')
    if value < input_type['min']:
        raise ntf.Notification.from_name('ValueTooSmall', input_type['min'])
    if value > input_type['max']:
        raise ntf.Notification.from_name('ValueTooBig', input_type['max'])
    return True


def check_double(value, input_type):
    """Checks if value is a real number within given range."""
    if not isinstance(value, (int, float)):
        raise ntf.Notification.from_name('ValidationTypeError', 'Double')
    if value < input_type['min']:
        raise ntf.Notification.from_name('ValueTooSmall', input_type['min'])
    if value > input_type['max']:
        raise ntf.Notification.from_name('ValueTooBig', input_type['max'])
    return True


def check_bool(value, input_type):
    """Checks if value is a boolean."""
    if not isinstance(value, bool):
        raise ntf.Notification.from_name('ValidationTypeError', 'Bool')
    return True


def check_string(value, input_type):
    """Checks if value is a string."""
    if not isinstance(value, str):
        raise ntf.Notification.from_name('ValidationTypeError', 'String')
    return True


def check_selection(value, input_type):
    """Checks if value is a valid option in given selection."""
    if value in input_type['values']:
        return True
    else:
        raise ntf.Notification.from_name('InvalidSelectionOption', input_type['name'], value)


def check_filename(value, input_type):
    """Placeholder for FileName validation."""
    return check_string(value, input_type)


def check_array(value, input_type):
    """Checks if value is an array of size within a given range."""
    if not isinstance(value, (list, str)):
        raise ntf.Notification.from_name('ValidationTypeError', 'Array')
    if len(value) < input_type['min']:
        raise ntf.Notification.from_name('NotEnoughItems', input_type['min'])
    elif len(value) > input_type['max']:
        raise ntf.Notification.from_name('TooManyItems', input_type['max'])
    return True


def check_record_key(children_keys, key, input_type):
    """Checks a single key within a record."""
    # if key is not found in specifications, it is considered to be valid
    if key not in input_type['keys']:
        raise ntf.Notification.from_name('UnknownRecordKey', key, input_type['type_name'])

    try:
        key_type = input_type['keys'][key]['default']['type']
    except KeyError:
        pass  # if default or type isn't specified, skip
    else:
        if key_type == 'obligatory':
            if key not in children_keys:
                raise ntf.Notification.from_name('MissingObligatoryKey', key,
                                                 input_type['type_name'])
    return True


def get_abstractrecord_type(node, input_type):
    """
    Returns the concrete TYPE of abstract record. ValidationErrors
    can occur if it is impossible to resolve the type.
    """
    try:
        type_node = node.type
    except AttributeError:
        raise ntf.Notification.from_name('ValidationTypeError', 'AbstractRecord')

    if type_node is None:
        try:
            concrete_type = input_type['default_descendant']
        except KeyError:
            raise ntf.Notification.from_name('MissingAbstractRecordType')
    else:
        try:
            concrete_type = input_type['implementations'][type_node.value]
        except KeyError:
            raise ntf.Notification.from_name('InvalidAbstractRecordType', type_node.value,
                                             input_type['name'])
    if concrete_type is None:
        raise ntf.Notification.from_name('MissingAbstractRecordType')
    return concrete_type
