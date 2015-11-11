"""
Customized pyyaml resolver
"""

import re


def resolve_scalar_tag(value):
    """Resolves a tag for scalar value."""
    if value == '':
        resolvers = _YAML_IMPLICIT_RESOLVERS.get('', [])
    else:
        resolvers = _YAML_IMPLICIT_RESOLVERS.get(value[0], [])
    resolvers += _YAML_IMPLICIT_RESOLVERS.get(None, [])
    for tag, regexp in resolvers:
        if regexp.match(value):
            return tag
    return _DEFAULT_SCALAR_TAG


_DEFAULT_SCALAR_TAG = 'tag:yaml.org,2002:str'

_YAML_IMPLICIT_RESOLVERS = {}


def _add_implicit_resolver(tag, regexp, first):
    """Adds resolver to implicit resolvers."""
    if first is None:
        first = [None]
    for char in first:
        _YAML_IMPLICIT_RESOLVERS.setdefault(char, []).append((tag, regexp))


_add_implicit_resolver(
    'tag:yaml.org,2002:bool',
    re.compile(r'''^(?:yes|Yes|YES|no|No|NO
                |true|True|TRUE|false|False|FALSE
                |on|On|ON|off|Off|OFF)$''', re.X),
    list('yYnNtTfFoO'))

_add_implicit_resolver(
    'tag:yaml.org,2002:int',
    re.compile(r'''^(?:[-+]?0b[0-1_]+
                |[-+]?0[0-7_]+
                |[-+]?(?:0|[1-9][0-9_]*)
                |[-+]?0x[0-9a-fA-F_]+
                |[-+]?[1-9][0-9_]*(?::[0-5]?[0-9])+)$''', re.X),
    list('-+0123456789'))

_add_implicit_resolver(
    'tag:yaml.org,2002:float',
    re.compile(r'''^(?:[-+]? ( \. [0-9]+ | [0-9]+ ( \. [0-9]* )? ) ( [eE] [-+]? [0-9]+ )?
                |\.[0-9_]+(?:[eE][-+][0-9]+)?
                |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\.[0-9_]*
                |[-+]?\.(?:inf|Inf|INF)
                |\.(?:nan|NaN|NAN))$''', re.X),
    list('-+0123456789.'))

_add_implicit_resolver(
    'tag:yaml.org,2002:merge',
    re.compile(r'^(?:<<)$'),
    ['<'])

_add_implicit_resolver(
    'tag:yaml.org,2002:null',
    re.compile(r'''^(?: ~
                |null|Null|NULL
                | )$''', re.X),
    ['~', 'n', 'N', ''])

_add_implicit_resolver(
    'tag:yaml.org,2002:timestamp',
    re.compile(r'''^(?:[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]
                |[0-9][0-9][0-9][0-9] -[0-9][0-9]? -[0-9][0-9]?
                 (?:[Tt]|[ \t]+)[0-9][0-9]?
                 :[0-9][0-9] :[0-9][0-9] (?:\.[0-9]*)?
                 (?:[ \t]*(?:Z|[-+][0-9][0-9]?(?::[0-9][0-9])?))?)$''', re.X),
    list('0123456789'))
