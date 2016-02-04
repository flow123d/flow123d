#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

from ist.globals import Globals


class InputType(object):
    SELECTION = 1
    RECORD = 2
    ABSTRACT = 4
    ABSTRACTRECORD = 4
    ARRAY = 8
    INTEGER = 16
    DOUBLE = 32
    PARAMETER = 64
    STRING = 128
    BOOL = 256
    FILENAME = 512

    MAIN_TYPE = SELECTION | RECORD | ABSTRACT

    def __eq__(self, other):
        if type(other) is int:
            return bool(self.value & other)

        if type(other) is str:
            return bool(self.value & getattr(self, other.upper(), 0))

        raise Exception('Unsupported comparison with type %s' % str(type(other)))

    def __init__(self):
        self.value = 0

    def __repr__(self):
        for value in dir(self):
            if value == value.upper():
                if getattr(self, value, 0) == self.value:
                    return '{value}'.format(value=value)
        return 'UNKNOWN'

    def parse(self, value):
        self.value = getattr(self, value.upper())
        return self


class Field(object):
    def __init__(self, names, t=str, index=False, subtype=None, save_as=None, required=False):
        """
        :type subtype: class
        """
        self.names = names if type(names) is list else [names]
        self.type = t
        self.index = index
        self.required = required
        self.subtype = subtype
        self.save_as = save_as or self.names[0]

    def parse(self, json_data):
        for name in self.names:
            if name in json_data:
                # parsable classes will handle parsing by themselves
                if callable(getattr(self.type, 'parse', None)):
                    if callable(getattr(self.type, 'set_subtype', None)):
                        instance = self.type()
                        instance.set_subtype(self.subtype)
                        return instance.parse(json_data[name])
                    return self.type().parse(json_data[name])

                return json_data[name]

        if self.required:
            raise

    @staticmethod
    def factory(*args, **kwargs):
        def create():
            return Field(*args, **kwargs)

        return create


class Parsable(object):
    __fields__ = []

    def parse(self, json_data={ }):
        for field in self.__fields__:
            value = field.parse(json_data)

            if field.index:
                Globals.items[value] = self

            self.__setattr__(field.save_as, value)
        return self

    def __repr__(self):
        if getattr(self, 'id', None):
            if getattr(self, 'name', None) and getattr(self, 'input_type', None):
                return "<{self.input_type} {self.name} [{self.id}]>".format(self=self)
            return "<{self.id}>".format(self=self)
        return super(Parsable, self).__repr__()

    def include_in_format(self):
        input_type = getattr(self, 'input_type', None)
        return input_type is not None and input_type == InputType.MAIN_TYPE


    @property
    def debug_name(self):
        # if getattr(self, 'name', None) and getattr(self, 'id', None):
        #     return '{self.name}[{self.id}]'.format(self=self)
        return self.name


class List(list):
    def __init__(self):
        super(List, self).__init__()
        self.subtype = str

    def set_subtype(self, subtype):
        self.subtype = subtype

    def parse(self, json_data):
        if self.subtype:
            if callable(getattr(self.subtype, 'parse', None)):
                for item in json_data:
                    self.append(self.subtype().parse(item))
            else:
                for item in json_data:
                    self.append(self.subtype(item))
        else:
            for item in json_data:
                self.append(item)
        return self

    def is_single(self):
        return len(self) == 1


class Dict(dict):
    def parse(self, json_data):
        for key, value in json_data.items():
            self[key] = value
        return self

    def __getattr__(self, item):
        return self.get(item, None)

class NotImplementedException(Exception):
    pass
