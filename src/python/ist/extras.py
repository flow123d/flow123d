#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

from __future__ import absolute_import
from ist.globals import Globals, FormatMode
from ist.base import Parsable, Field, List, Dict, InputType
from ist.utils.htmltree import htmltree
from ist.utils.texlist2 import TexList


class TypeReference(Parsable):
    """
    :type reference  : unicode
    """
    __fields__ = [
    ]

    def parse(self, json_data={ }):
        self.reference = json_data
        return self

    def __init__(self):
        self.reference = None

    def get_reference(self):
        """
        :rtype : ist.nodes.TypeRecord or ist.nodes.TypeSelection or ist.nodes.TypeAbstract
        or ist.nodes.TypeString or nodes.Integer or ist.nodes.TypeDouble or ist.nodes.TypeArray
        or ist.nodes.TypeParameter or ist.nodes.TypeFilename or ist.nodes.TypeBool
        """
        return Globals.items[self.reference]

    @property
    def target(self):
        return self.get_reference()

    def __repr__(self):
        ref = self.get_reference().input_type
        return '{self.reference}={ref}'.format(self=self, ref=ref)


class TypeSelectionValue(Parsable):
    """
    :type name           : unicode
    :type description    : unicode
    :type parent         : Parsable
    """
    __fields__ = [
        Field('name'),
        Field('description'),
    ]

    def __init__(self):
        super(TypeSelectionValue, self).__init__()
        self.name = None
        self.description = None

    @property
    def href_id(self):
        if self.parent:
            parent_id = str(self.parent.href_id)
            if FormatMode.format_mode == FormatMode.LATEX_MODE:
                if parent_id.startswith('IT::'):
                    return '{}::{}'.format(parent_id[4:], TexList.name_mode(self.name))
                return '{}::{}'.format(parent_id, TexList.name_mode(self.name))

            return '{self.parent.href_id}-{self.name}'.format(self=self)
        return self.name


    @property
    def href_name(self):
        return self.name


class TypeRecordKeyDefault(Parsable):
    """
    :type type           : unicode
    :type value          : unicode
    """
    __fields__ = [
        Field('type'),
        Field('value'),
    ]

    def __init__(self):
        super(TypeRecordKeyDefault, self).__init__()
        self.type = None
        self.value = None


class TypeRecordKey(Parsable):
    """
    :type key            : unicode
    :type type           : ist.extras.TypeReference
    :type default        : ist.extras.TypeRecordKeyDefault
    :type description    : unicode
    """
    __fields__ = [
        Field('key'),
        Field('type', t=TypeReference),
        Field('default', t=TypeRecordKeyDefault, link_to_parent=True),
        Field('description'),
    ]

    def __init__(self):
        super(TypeRecordKey, self).__init__()
        self.key = None
        self.type = None
        self.default = None
        self.description = None

    def include_in_format(self):
        return True

    @property
    def href_id(self):
        if self.parent:
            parent_id = str(self.parent.href_id)
            if FormatMode.format_mode == FormatMode.LATEX_MODE:
                if parent_id.startswith('IT::'):
                    return '{}::{}'.format(parent_id[4:], TexList.name_mode(self.key))
                return '{}::{}'.format(parent_id, TexList.name_mode(self.key))

            return '{self.parent.href_id}-{self.key}'.format(self=self)
        return self.key

    @property
    def href_name(self):
        return htmltree.secure(self.key)


class TypeRange(Parsable):
    __fields__ = []

    replacements = {
        '2147483647': 'INT32 MAX',
        '4294967295': 'UINT32 MAX',
        '-2147483647': 'INT32 MIN',
        '1.79769e+308': '+inf',
        '-1.79769e+308': '-inf',
        '': 'unknown range'
    }

    def parse(self, json_data={ }):
        self.min = json_data[0]
        self.max = json_data[1]
        return self

    def __init__(self):
        self.min = ''
        self.max = ''
        self.always_visible = True

    def is_pointless(self):
        """
         Whether is information within this instance beneficial
        """
        return self._format() in ('[0, ]', '[, ]', '(-inf, +inf)')

    def _format(self):
        """
        Method will will return string representation of this range
        :return:
        """
        min_value = self.replacements.get(str(self.min), str(self.min))
        max_value = self.replacements.get(str(self.max), str(self.max))
        l_brace = '(' if min_value.find('inf') != -1 else '['
        r_brace = ')' if max_value.find('inf') != -1 else ']'

        return '{l_brace}{min_value}, {max_value}{r_brace}'.format(
            l_brace=l_brace, r_brace=r_brace,
            min_value=min_value, max_value=max_value)

    def __repr__(self):
        """
        method will return string representation if is meaningful or flag always visible
        is True
        """
        return self._format() if self.always_visible or not self.is_pointless() else ''


class TypeAttributeParameter(Parsable):
    """
    :type name           : unicode
    :type reference      : ist.extras.TypeReference
    """
    __fields__ = [
    ]

    def __init__(self):
        self.name = None
        self.reference = None

    def parse(self, json_data={ }):
        item = json_data.items()[0]
        self.name = str(item[0])
        self.reference = TypeReference().parse(item[1])
        return self

    def __repr__(self):
        if self.name and self.reference:
            return '<{self.name} -> {self.reference}>'.format(self=self)
        return '<>'


class TypeAttributes(Parsable):
    """
    :type obsolete       : unicode
    :type link_name      : unicode
    :type parameters     : list[TypeAttributeParameter]
    :type generic_type   : ist.extras.TypeReference
    :type input_type     : InputType
    """
    __fields__ = [
        Field('obsolete', t=str),
        Field('link_name', index=True),
        Field('parameters', t=List, subtype=TypeAttributeParameter),
        Field('generic_type', t=TypeReference),
    ]

    def __init__(self):
        super(TypeAttributes, self).__init__()
        self.obsolete = None
        self.link_name = None
        self.parameters = None
        self.generic_type = None
        self.input_type = InputType().parse('')

    def __repr__(self):
        if self.obsolete is None and self.link_name is None \
                and self.parameters is None and self.generic_type is None:
            return '{}'
        return super(TypeAttributes, self).__repr__()

    def get_parameters_dict(self):
        """
        :rtype : Dict
        """
        if not self.parameters:
            return Dict()

        result = Dict()
        for p in self.parameters:
            result[p.name] = p
        return result