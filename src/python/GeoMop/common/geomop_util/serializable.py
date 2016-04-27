"""Serialization module.

.. codeauthor:: Tomas Krizek <tomas.krizek1@tul.cz>
"""

import copy
from enum import Enum


class Serializable:
    """Class defines special operations during the serialization process.

    It can:
      - exclude certain keys from serialization
      - delete key from serialized file (same as excluded, but used for deprecated / removed keys)
      - set default values for keys if they are not specified in the source data
      - allow serialization of nested objects

    Serialization of nested objects.
    Define nested keys in the composite dictionary. As a value, pass in the class to be instanced.
    If the class needs to reference itself, you can define __serializable__ after the class
    definition. See testing/common/test_serializable.py for example.
    """

    def __init__(self, excluded=None, deleted=None, default=None, composite=None):
        self.excluded = excluded if excluded is not None else []
        self.deleted = deleted if deleted is not None else []
        self.default = default if default is not None else {}
        self.composite = composite if composite is not None else {}

    @staticmethod
    def load(data, cls=None):
        """Create object data structure from native dict."""
        if cls is not None:
            if hasattr(cls, '__serializable__'):
                rules = cls.__serializable__
            else:
                rules = Serializable()
        else:
            # nothing to do, no rules defined
            return data

        if isinstance(data, list):
            deserialized = []
            for item in data:
                deserialized.append(Serializable.load(item, cls))
            return deserialized

        if data is None:
            return cls()
        elif not isinstance(data, dict):
            return cls(data)

        for exclude in (rules.excluded + rules.deleted):
            if exclude in data:
                del data[exclude]

        for key, value in rules.default.items():
            if key not in data:
                data[key] = value

        # __all__: set default composite
        composite = {}
        if '__all__' in rules.composite:
            default_type = rules.composite['__all__']
            composite = {key: default_type for key in data}
        # override default composite
        composite.update(rules.composite)

        # recursively resolve composite
        for key, class_ in composite.items():
            if key in data:
                subdata = data[key]
                data[key] = Serializable.load(subdata, class_)

        # finally, construct the class
        return cls(**data)

    @staticmethod
    def dump(data):
        """Create serializable data structure from provided data."""
        if hasattr(data, '__serializable__'):
            rules = data.__serializable__
        elif hasattr(data, '__dict__'):
            rules = Serializable()
        elif isinstance(data, list):
            serialized = []
            for item in data:
                serialized.append(Serializable.dump(item))
            return serialized
        else:
            # nothing to do, no rules defined
            return data

        # different serialization for dict, enum and class
        if isinstance(data, dict):
            out = dict(copy.copy(data))
        elif isinstance(data, Enum):
            return data.value
        else:
            out = copy.copy(data.__dict__)

        for exclude in (rules.excluded + rules.deleted):
            if exclude in out:
                del out[exclude]

        # __all__: set default composite
        composite = {}
        if '__all__' in rules.composite:
            default_type = rules.composite['__all__']
            composite = {key: default_type for key in out}
        # override default composite
        composite.update(rules.composite)

        # recursively resolve composite
        for key, class_ in composite.items():
            if key in out:
                subdata = out[key]
                out[key] = Serializable.dump(subdata)

        return out
