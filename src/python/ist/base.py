#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

from ist.globals import Globals, FormatMode
from ist.utils.htmltree import htmltree
from ist.utils.texlist2 import TexList
from utils.logger import Logger
import collections


class InputType(object):
    """
    Class InputType is enum-like class for types
    """

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
    TUPLE = 1024

    MAIN_TYPE = SELECTION | RECORD | ABSTRACT | TUPLE

    def __eq__(self, other):
        if other is None:
            return False

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
    """
    Class Field registers new field on class
    """

    def __init__(self, names, t=str, index=False, subtype=None, save_as=None, required=False, link_to_parent=False, default=None):
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
        self.default = default

    def parse(self, json_data):
        for name in self.names:
            if name in json_data:
                # parsable classes will handle parsing by themselves
                if isinstance(getattr(self.type, 'parse', None), collections.abc.Callable):
                    instance = self.type()
                    if isinstance(getattr(self.type, 'set_subtype', None), collections.abc.Callable):
                        instance.set_subtype(self.subtype)
                    return instance.parse(json_data[name])

                return json_data[name]

        if self.default:
            return self.default

        if self.required:
            raise

    @staticmethod
    def factory(*args, **kwargs):
        def create():
            return Field(*args, **kwargs)

        return create


class Parsable(object):
    """
    Class Parsable is abstract helper class for parsable objects
    """

    __fields__ = []

    def __init__(self):
        self.parent = None
        self.unique_name = None
        self.references = list()
        self.secnerefer = list()

    def add_link(self, target):
        """
        :type target: Parsable
        """
        left = self.get_generic_root()
        right = target.get_generic_root()
        left.add_ref(right)
        right.add_fer(left)

    def add_ref(self, ref):
        if ref.input_type == InputType.MAIN_TYPE:
            self.references.append(ref)

    def add_fer(self, fer):
        if fer.input_type == InputType.MAIN_TYPE:
            self.secnerefer.append(fer)

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
                    if (self.has_generic_link()):
                        value="instance_"+value
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
        name = getattr(self, 'name', '')

        if not (input_type == InputType.MAIN_TYPE):
            Logger.instance().info('[SKIP] not main type %s' % self)
            return False

        if self.has_generic_link():
            Logger.instance().info('[SKIP] has generic link \n%s \npoints to \n    %s' % (self, self.get_generic_root()))
            return False
        
        if name == 'EmptyRecord':
            Logger.instance().info('[SKIP] name is EmptyRecord %s' % self)
            return False
        return True

    def has_generic_link(self):
        return getattr(self, 'generic_type', None) is not None

    def get_parameter_dict(self):
        return getattr(self, 'parameters', {})

    def get_generic_root(self):
        """
        return generic root of this object
        :rtype: Parsable
        """
        root = self
        while True:
            try:
                if root.generic_type:
                    root = root.generic_type.get_reference()
                else:
                    break
            except: break
        return root


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

        # return self.href_id
        return getattr(self, 'name', None)

    @property
    def href_id(self):
        """
        Method will return unique id of this item
        if link_name was specified in attributes, it will be used
        :return:
        """
        if getattr(self, 'attributes', None) and self.attributes.link_name:
            value = self.attributes.link_name
        else:
            value = self.get('unique_name', 'name', 'id')

        if FormatMode.format_mode == FormatMode.LATEX_MODE:
            return 'IT::{}'.format(TexList.name_mode(value))
        return htmltree.secure(value)

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
    """
    Class List adds extra functionality to classic list (such as parsing)
    """

    def __init__(self):
        super(List, self).__init__()
        self.subtype = str

    def set_subtype(self, subtype):
        self.subtype = subtype

    def parse(self, json_data):
        if self.subtype:
            if isinstance(getattr(self.subtype, 'parse', None), collections.abc.Callable):
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


class SmartList(List):
    """
    Class SmartList adds extra functionality to List
    parsing is compatible with previous versions
    """

    def parse(self, json_data):
        if type(json_data) is dict:
            return super(SmartList, self).parse([dict([x]) for x in json_data.items()])
        return super(SmartList, self).parse(json_data)


class Dict(dict):
    """
    Class Dict adds extra functionality to classic dict
    """

    def parse(self, json_data):
        for key, value in json_data.items():
            self[key] = value
        return self

    def __getattr__(self, item):
        return self.get(item, None)


class Unicode(Parsable):
    """
    Class Unicode represents str and unicode data
    :type value          : unicode
    """

    def __init__(self):
        super(Parsable, self).__init__()
        self.value = ''

    def parse(self, json_data=''):
        self.value = str(json_data)
        return self

    def __str__(self):
        return str.__str__(str(self.value))

    def __repr__(self):
        return str.__repr__(str(self.value))

    def __bool__(self):
        return bool(self.value)

    @property
    def href_name(self):
        # if self.parent:
        # return '{self.parent.href_name}-{self.value}'.format(self=self)
        return self.value

    @property
    def href_id(self):
        if self.parent:
            parent_id = str(self.parent.href_id)
            if FormatMode.format_mode == FormatMode.LATEX_MODE:
                if parent_id.startswith('IT::'):
                    return '{}::{}'.format(parent_id[4:], TexList.name_mode(self.value))
                return '{}::{}'.format(parent_id, TexList.name_mode(self.value))

            return '{self.parent.href_id}:{self.value}'.format(self=self)
        return self.value


class NotImplementedException(Exception):
    """
    Simple exception class
    """
    pass
