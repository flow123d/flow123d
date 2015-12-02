"""
Customized PyYAML Constructor
"""

import re
import datetime


def construct_scalar(value, tag):
    """Constructs a scalar value of the correct python type."""
    if tag not in _CONSTRUCTORS:
        raise Exception('Tag {tag} is not supported for scalar values.'
                        .format(tag=tag))
    return _CONSTRUCTORS[tag](value)


# pylint: disable=unused-argument
def _construct_yaml_null(value):
    """YAML constructor for null."""
    return None


_BOOL_VALUES = {
    'yes':      True,
    'no':       False,
    'true':     True,
    'false':    False,
    'on':       True,
    'off':      False,
}


def _construct_yaml_bool(value):
    """YAML constructor for bool."""
    return _BOOL_VALUES[value.lower()]


def _construct_yaml_int(value):
    """YAML constructor for integer."""
    sign = +1
    if value[0] == '-':
        sign = -1
    if value[0] in '+-':
        value = value[1:]
    if value == '0':
        return 0
    elif value.startswith('0b'):
        return sign*int(value[2:], 2)
    elif value.startswith('0x'):
        return sign*int(value[2:], 16)
    elif value[0] == '0':
        return sign*int(value, 8)
    elif ':' in value:
        digits = [int(part) for part in value.split(':')]
        digits.reverse()
        base = 1
        value = 0
        for digit in digits:
            value += digit*base
            base *= 60
        return sign*value
    else:
        return sign*int(value)


def _construct_yaml_float(value):
    """YAML constructor for float."""
    value = value.lower()
    sign = +1
    if value[0] == '-':
        sign = -1
    if value[0] in '+-':
        value = value[1:]
    if value == '.inf':
        return sign*float('inf')
    elif value == '.nan':
        return float('NaN')
    elif ':' in value:
        digits = [float(part) for part in value.split(':')]
        digits.reverse()
        base = 1
        value = 0.0
        for digit in digits:
            value += digit*base
            base *= 60
        return sign*value
    else:
        return sign*float(value)


_TIMESTAMP_REGEX = re.compile((
    r'^(?P<year>[0-9][0-9][0-9][0-9])'
    r'-(?P<month>[0-9][0-9]?)'
    r'-(?P<day>[0-9][0-9]?)'
    r'-(?:(?:[Tt]|[ \t]+)'
    r'-(?P<hour>[0-9][0-9]?)'
    r'-:(?P<minute>[0-9][0-9])'
    r'-:(?P<second>[0-9][0-9])'
    r'-(?:\.(?P<fraction>[0-9]*))?'
    r'-(?:[ \t]*(?P<tz>Z|(?P<tz_sign>[-+])(?P<tz_hour>[0-9][0-9]?)'
    r'-(?::(?P<tz_minute>[0-9][0-9]))?))?)?$'), re.X)


def _construct_yaml_timestamp(value):
    """YAML constructor for timestamp."""
    match = _TIMESTAMP_REGEX.match(value)
    values = match.groupdict()
    year = int(values['year'])
    month = int(values['month'])
    day = int(values['day'])
    if not values['hour']:
        return datetime.date(year, month, day)
    hour = int(values['hour'])
    minute = int(values['minute'])
    second = int(values['second'])
    fraction = 0
    if values['fraction']:
        fraction = values['fraction'][:6]
        while len(fraction) < 6:
            fraction += '0'
        fraction = int(fraction)
    delta = None
    if values['tz_sign']:
        tz_hour = int(values['tz_hour'])
        tz_minute = int(values['tz_minute'] or 0)
        delta = datetime.timedelta(hours=tz_hour, minutes=tz_minute)
        if values['tz_sign'] == '-':
            delta = -delta
    data = datetime.datetime(year, month, day, hour, minute, second, fraction)
    if delta:
        data -= delta
    return data


def _construct_yaml_str(value):
    """YAML constructor for string."""
    return str(value)


_CONSTRUCTORS = {
    'tag:yaml.org,2002:null': _construct_yaml_null,
    'tag:yaml.org,2002:bool': _construct_yaml_bool,
    'tag:yaml.org,2002:int': _construct_yaml_int,
    'tag:yaml.org,2002:float': _construct_yaml_float,
    'tag:yaml.org,2002:timestamp': _construct_yaml_timestamp,
    'tag:yaml.org,2002:str': _construct_yaml_str
}

