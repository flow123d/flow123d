#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

from __future__ import absolute_import
from ist.globals import Globals
from ist.utils.htmltree import htmltree
from utils.logger import Logger


class InputType(object):
    UNKNOWN = 0
    SELECTION = 4
    RECORD = 1
    ABSTRACT = 2
    ABSTRACT_RECORD = ABSTRACT
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
        self.value = getattr(self, value.upper(), InputType.UNKNOWN)
        return self


class Field(object):
    def __init__(self, names, t=str, index=False, subtype=None, save_as=None, required=False, link_to_parent=False):
        """
        :type subtype: class
        """
        self.names = names if type(names) is list else [names]
        self.type = t
        self.index = index
        self.required = required
        self.subtype = subtype
        self.save_as = save_as or self.names[0]
        self.link_to_parent = link_to_parent

    def parse(self, json_data):
        for name in self.names:
            if name in json_data:
                # parsable classes will handle parsing by themselves
                if callable(getattr(self.type, 'parse', None)):
                    instance = self.type()
                    if callable(getattr(self.type, 'set_subtype', None)):
                        instance.set_subtype(self.subtype)
                    return instance.parse(json_data[name])

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

    def __init__(self):
        self.parent = None
        self.unique_name = None
        self.references = list()

    def add_ref(self, ref):
        self.references.append(ref)

    def parse(self, json_data={}):
        for field in self.__fields__:
            value = field.parse(json_data)

            if field.link_to_parent:
                if getattr(value, '__iter__', None):
                    for v in value:
                        if v:
                            v.set_parent(self)
                else:
                    if value:
                        value.set_parent(self)

            if field.index:
                if value:
                    unique_name = Globals.save(value, self)
                    if value != unique_name:
                        self.unique_name = unique_name

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
        # return '{self.name}[{self.id}]'.format(self=self)
        return self.name

    @property
    def href_name(self):
        """
        Method will return unique name of this item
        if link_name was specified in attributes, it will be used
        :return:
        """
        # if getattr(self, 'attributes', None) and self.attributes.link_name:
        # return self.attributes.link_name

        # if getattr(self, 'unique_name', None):
        #     return self.unique_name

        if getattr(self, 'name', None):
            return self.name


    @property
    def href_id(self):
        """
        Method will return unique id of this item
        if link_name was specified in attributes, it will be used
        :return:
        """
        # if getattr(self, 'attributes', None) and self.attributes.link_name:
        # return self.attributes.link_name
        #
        if getattr(self, 'unique_name', None):
            return htmltree.secure(self.unique_name)
        #
        if getattr(self, 'name', None):
            return htmltree.secure(self.name)

        if getattr(self, 'id', None):
            return htmltree.secure(self.id)

    def get(self, *args):
        """
        Method will return first matching valid property of this object
        :param args: list[str]
        :return:
        """
        for arg in args:
            value = getattr(self, arg, None)
            if value:
                return value
        Logger.instance().debug('Cannot find {args} on {self}'.format(args=args, self=self))
        return None

    def gets(self, *args):
        """
        Method will return list of values by given prop names
        :param args:
        :return: list
        """
        return [self.get(arg) for arg in args]

    def get_fields(self, *args):
        """
        Method will return first matching property from its sub keys
        :param args: list[str]
        :return:
        """
        raise NotImplementedException('get_fields invoked on simple object')

    def set_parent(self, parent):
        """
        :type parent: Parsable
        """
        self.parent = parent

    def get_references(self):
        return list(set(self.references))
        # result = []
        # for item in Globals.iterate():
        # for key in getattr(item, 'keys', []):
        #         if key.type.get_reference().id == self.id:
        #             result.append(item)
        #         if key.type.get_reference().input_type == InputType.ARRAY:
        #             if key.type.get_reference().subtype.get_reference().id == self.id:
        #                 result.append(item)
        #
        #     for imp in getattr(item, 'implementations', []):
        #         if imp.get_reference().id == self.id:
        #             result.append(item)
        #
        # return list(set(result))


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


class Unicode(Parsable):
    """
    :type value          : unicode
    """

    def __init__(self):
        super(Parsable, self).__init__()
        self.value = u''

    def parse(self, json_data=u''):
        self.value = unicode(json_data)
        return self

    def __str__(self):
        return unicode.__str__(unicode(self.value))

    def __repr__(self):
        return unicode.__repr__(unicode(self.value))

    def __nonzero__(self):
        return bool(self.value)

    @property
    def href_name(self):
        # if self.parent:
        # return '{self.parent.href_name}-{self.value}'.format(self=self)
        return self.value

    @property
    def href_id(self):
        if self.parent:
            return '{self.parent.href_id}-{self.value}'.format(self=self)
        return self.value


class NotImplementedException(Exception):
    pass
