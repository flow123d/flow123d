#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

from ist.extras import TypeSelectionValue, TypeReference, TypeRecordKey, TypeRange, TypeAttributes
from ist.base import Field, Parsable, Dict, InputType, List, Unicode


class TypeSelection(Parsable):
    """
    Class defining "Selection" type in IST
    :type id                 : unicode
    :type values             : list[ist.extras.TypeSelectionValue]
    :type name               : unicode
    :type input_type         : InputType
    :type attributes         : ist.extras.TypeAttributes
    :type description        : unicode
    """
    __fields__ = [
        Field("id", index=True),
        Field("values", t=List, subtype=TypeSelectionValue, link_to_parent=True),
        Field(["name", "type_name"], index=True),
        Field("input_type", t=InputType),
        Field("attributes", t=TypeAttributes),
        Field("description"),
    ]

    def __init__(self):
        super(TypeSelection, self).__init__()
        self.id = None
        self.values = None
        self.name = None
        self.input_type = None
        self.attributes = None
        self.description = None

    def include_in_format(self):
        return self.name.find('TYPE') == -1

    def get_fields(self, *args):
        for arg in args:
            for sub_item in self.values:
                if sub_item.name.lower() == arg.lower():
                    return sub_item


class TypeRecord(Parsable):
    """
    Class defining "Record" type in IST
    :type id                 : unicode
    :type keys               : list[ist.extras.TypeRecordKey]
    :type name               : unicode
    :type implements         : List[ist.extras.TypeReference]
    :type input_type         : InputType
    :type attributes         : ist.extras.TypeAttributes
    :type description        : unicode
    :type reducible_to_key   : Unicode
    """
    __fields__ = [
        Field("id", index=True),
        Field("keys", t=List, subtype=TypeRecordKey, link_to_parent=True),
        Field(["name", "type_name"], index=True),
        Field("implements", t=List, subtype=TypeReference),
        Field("input_type", t=InputType),
        Field("attributes", t=TypeAttributes),
        Field("description"),
        Field("reducible_to_key", t=Unicode, link_to_parent=True),
    ]

    def __init__(self):
        super(TypeRecord, self).__init__()
        self.id = None
        self.keys = None
        self.name = None
        self.implements = None
        self.input_type = None
        self.attributes = None
        self.description = None
        self.reducible_to_key = None

    def get_fields(self, *args):
        for arg in args:
            for sub_item in self.keys:
                if sub_item.key.lower() == arg.lower():
                    return sub_item


class TypeAbstract(Parsable):
    """
    Class defining "Abstract" type in IST
    :type id                 : unicode
    :type name               : unicode
    :type input_type         : InputType
    :type attributes         : ist.extras.TypeAttributes
    :type description        : unicode
    :type implementations    : List[ist.extras.TypeReference]
    :type default_descendant : unicode
    """
    __fields__ = [
        Field("id", index=True),
        Field(["name", "type_name"], index=True),
        Field("input_type", t=InputType),
        Field("attributes", t=TypeAttributes),
        Field("description"),
        Field("implementations", t=List, subtype=TypeReference),
        Field("default_descendant", t=TypeReference),
    ]

    def __init__(self):
        super(TypeAbstract, self).__init__()
        self.id = None
        self.name = None
        self.input_type = None
        self.attributes = None
        self.description = None
        self.implementations = None
        self.default_descendant = None


class TypeString(Parsable):
    """
    Class defining "String" type in IST
    :type id                 : unicode
    :type name               : unicode
    :type input_type         : InputType
    :type attributes         : ist.extras.TypeAttributes
    """
    __fields__ = [
        Field("id", index=True),
        Field(["name", "type_name"]),
        Field("input_type", t=InputType),
        Field("attributes", t=TypeAttributes),
    ]

    def __init__(self):
        super(TypeString, self).__init__()
        self.id = None
        self.name = None
        self.input_type = None
        self.attributes = None


class TypeDouble(Parsable):
    """
    Class defining "Double" type in IST
    :type id                 : unicode
    :type range              : ist.extras.TypeRange
    :type name               : unicode
    :type input_type         : InputType
    :type attributes         : ist.extras.TypeAttributes
    """
    __fields__ = [
        Field("id", index=True),
        Field("range", t=TypeRange),
        Field(["name", "type_name"]),
        Field("input_type", t=InputType),
        Field("attributes", t=TypeAttributes),
    ]

    def __init__(self):
        super(TypeDouble, self).__init__()
        self.id = None
        self.range = None
        self.name = None
        self.input_type = None
        self.attributes = None


class TypeFilename(Parsable):
    """
    Class defining "FileName" type in IST
    :type id                 : unicode
    :type name               : unicode
    :type file_mode          : unicode
    :type attributes         : ist.extras.TypeAttributes
    :type input_type         : InputType
    """
    __fields__ = [
        Field("id", index=True),
        Field(["name", "type_name"]),
        Field("file_mode"),
        Field("attributes", t=TypeAttributes),
        Field("input_type", t=InputType),
    ]

    def __init__(self):
        super(TypeFilename, self).__init__()
        self.id = None
        self.name = None
        self.file_mode = None
        self.attributes = None
        self.input_type = None


class TypeBool(Parsable):
    """
    Class defining "Bool" type in IST
    :type id                 : unicode
    :type name               : unicode
    :type input_type         : InputType
    :type attributes         : ist.extras.TypeAttributes
    """
    __fields__ = [
        Field("id", index=True),
        Field(["name", "type_name"]),
        Field("input_type", t=InputType),
        Field("attributes", t=TypeAttributes),
    ]

    def __init__(self):
        super(TypeBool, self).__init__()
        self.id = None
        self.name = None
        self.input_type = None
        self.attributes = None


class TypeInteger(Parsable):
    """
    Class defining "Integer" type in IST
    :type id                 : unicode
    :type range              : ist.extras.TypeRange
    :type name               : unicode
    :type input_type         : InputType
    :type attributes         : ist.extras.TypeAttributes
    """
    __fields__ = [
        Field("id", index=True),
        Field("range", t=TypeRange),
        Field(["name", "type_name"]),
        Field("input_type", t=InputType),
        Field("attributes", t=TypeAttributes),
    ]

    def __init__(self):
        super(TypeInteger, self).__init__()
        self.id = None
        self.range = None
        self.name = None
        self.input_type = None
        self.attributes = None


class TypeArray(Parsable):
    """
    Class defining "Array" type in IST
    :type id                 : unicode
    :type range              : ist.extras.TypeRange
    :type subtype            : ist.extras.TypeReference
    :type name               : unicode
    :type input_type         : InputType
    :type attributes         : ist.extras.TypeAttributes
    """
    __fields__ = [
        Field("id", index=True),
        Field("range", t=TypeRange),
        Field("subtype", t=TypeReference),
        Field(["name", "type_name"]),
        Field("input_type", t=InputType),
        Field("attributes", t=TypeAttributes),
    ]

    def __init__(self):
        super(TypeArray, self).__init__()
        self.id = None
        self.range = None
        self.subtype = None
        self.name = None
        self.input_type = None
        self.attributes = None


class TypeParameter(Parsable):
    """
    Class defining "Parameter" type in IST
    :type id                 : unicode
    :type name               : unicode
    :type input_type         : InputType
    :type attributes         : ist.extras.TypeAttributes
    """
    __fields__ = [
        Field("id", index=True),
        Field(["name", "type_name"]),
        Field("input_type", t=InputType),
        Field("attributes", t=TypeAttributes),
    ]

    def __init__(self):
        super(TypeParameter, self).__init__()
        self.id = None
        self.name = None
        self.input_type = None
        self.attributes = None
