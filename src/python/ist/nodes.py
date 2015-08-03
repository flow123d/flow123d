# encoding: utf-8
# author:   Jan Hybs
from ist.globals import Globals
from ist.utils.utils import Field, TypedList, AttributeDict


class ISTNode(object):
    # type of the object (similar to __class__.__name__)
    __type__ = ''

    # name of the field, which stores value by which can be this object searched
    __name_field__ = ''

    # name of the field, which stores value, which is the most important to this object
    __value_field__ = ''

    _fields = []


    def get_name(self):
        return getattr(self, self.__name_field__)


    def get_type(self):
        return self.__type__


    def get_value(self):
        return getattr(self, self.__value_field__)

    def parse(self, o={ }):
        for field in self._fields:
            # parse list
            if issubclass(field.type, TypedList) or issubclass(field.type, AttributeDict):
                value = field.value.copy().parse(o.get(field.name, []))
            # parse more complex objects
            elif issubclass(field.type, ISTNode):
                value = field.value.copy().parse(o.get(field.name, { }))
            # simple values
            else:
                value = o.get(field.name)

            # create attribute on class instance
            self.__setattr__(field.name, value)
            if field.name == 'id' and value is not None:
                Globals.items[value] = self

        return self

    def include_in_format(self):
        return True

    def get(self, *args):
        for arg in args:
            value = getattr(self, arg, None)
            if value is not None:
                return value

        raise Exception('no valid attribute within {} found on {}'.format(args, self.__class__.__name__))


    def copy(self):
        """
        Return copy of this instance
        :return:
        """
        return self.__class__()

    def __repr__(self):
        return '<{self.__class__.__name__} {self._fields}>'.format(self=self)


class Reference(ISTNode):
    __type__ = 'Reference'
    __value_field__ = 'ref_id'

    def parse(self, o=''):
        self.ref_id = o
        return self

    def get_reference(self):
        return Globals.items.get(str(self.ref_id), None)

    def copy(self):
        return Reference()

    def __repr__(self):
        return "<Reference ({self.ref_id}) -> {ref}>".format(self=self, ref=self.get_reference())

    def __nonzero__(self):
        return self.get_reference() is not None


class NumberRange(ISTNode):
    __type__ = 'Range'
    __name_field__ = ''
    __value_field__ = ''

    def __init__(self, always_visible=True):
        super(NumberRange, self).__init__()
        self.min = self.max = ''
        self.always_visible = always_visible

    replacements = {
        '2147483647': 'INT32 MAX',
        '4294967295': 'UINT32 MAX',
        '-2147483647': 'INT32 MIN',
        '1.79769e+308': '+inf',
        '-1.79769e+308': '-inf',
        '': 'unknown range'
    }

    def parse(self, o=[]):
        self.min = o[0] if len(o) > 0 else ''
        self.max = o[1] if len(o) > 1 else ''

        return self

    def is_pointless(self):
        """
        Wheather is information within this instance beneficial
        """
        return self._format() in ('[0, ]', '[, ]', '(-inf, +inf)')

    def _format(self):
        """
        Method will will return string representation of this instance
        :return:
        """
        min_value = self.replacements.get(str(self.min), str(self.min))
        max_value = self.replacements.get(str(self.max), str(self.max))
        l_brace = '(' if min_value.find('inf') != -1 else '['
        r_brace = ')' if max_value.find('inf') != -1 else ']'

        return '{l_brace}{min_value}, {max_value}{r_brace}'.format(
            l_brace=l_brace, r_brace=r_brace,
            min_value=min_value, max_value=max_value)

    def copy(self):
        return NumberRange(self.always_visible)

    def __repr__(self):
        """
        method will return string representation if is meaningful or flag always visible
        is True
        :return:
        """
        return self._format() if self.always_visible or not self.is_pointless() else ''


class DoubleRange(NumberRange):
    __type__ = 'Range'
    __name_field__ = ''

    def __init__(self, always_visible=False):
        super(DoubleRange, self).__init__(always_visible)

    def is_pointless(self):
        return self._format() == '[, ]'

    def copy(self):
        return DoubleRange(self.always_visible)


class AbstractNode(ISTNode):
    _fields = ISTNode._fields + [
        Field('id'),
        Field('input_type')
    ]


class Attribute(AbstractNode):
    __type__ = 'Attribute'
    __name_field__ = 'name'
    __value_field__ = 'value'

    _fields = AbstractNode._fields + [
        Field('name'),
        Field('value')
    ]


class DescriptionNode(AbstractNode):
    _fields = AbstractNode._fields + [
        Field('description')
    ]


class ComplexNode(DescriptionNode):
    _fields = DescriptionNode._fields + [
        Field('attributes', AttributeDict(Attribute))
    ]


class SelectionValue(DescriptionNode):
    __type__ = 'Selection value'
    __name_field__ = 'name'
    __value_field__ = 'name'

    _fields = DescriptionNode._fields + [
        Field('name')
    ]


class RecordKeyDefault(AbstractNode):
    __type__ = 'Defaults'
    __name_field__ = ''
    __value_field__ = 'value'

    _fields = AbstractNode._fields + [
        Field('type'),
        Field('value')
    ]


class RecordKey(DescriptionNode):
    __type__ = 'Record key'
    __name_field__ = 'key'
    __value_field__ = 'default'

    _fields = DescriptionNode._fields + [
        Field('key'),
        Field('type', Reference()),
        Field('default', RecordKeyDefault())
    ]

    def include_in_format(self):
        return self.key.find('TYPE') == -1


class String(AbstractNode):
    __type__ = 'String'
    __name_field__ = ''
    __value_field__ = 'name'

    _fields = AbstractNode._fields + [
        Field('name'),
        Field('full_name')
    ]


class Double(AbstractNode):
    __type__ = 'Double'
    __name_field__ = ''
    __value_field__ = 'range'

    _fields = AbstractNode._fields + [
        Field('name'),
        Field('full_name'),
        Field('range', DoubleRange())
    ]


class Integer(AbstractNode):
    __type__ = 'Integer'
    __name_field__ = ''
    __value_field__ = 'range'

    _fields = AbstractNode._fields + [
        Field('name'),
        Field('full_name'),
        Field('range', NumberRange())
    ]


class FileName(AbstractNode):
    __type__ = 'FileName'
    __name_field__ = ''
    __value_field__ = 'file_mode'

    _fields = AbstractNode._fields + [
        Field('name'),
        Field('full_name'),
        Field('file_mode')
    ]


class Bool(AbstractNode):
    __type__ = 'Bool'
    __name_field__ = ''
    __value_field__ = ''

    _fields = AbstractNode._fields + [
        Field('name'),
        Field('full_name')
    ]


class Array(AbstractNode):
    __type__ = 'Array'
    __name_field__ = ''
    __value_field__ = 'range'

    _fields = AbstractNode._fields + [
        Field('range', NumberRange(False)),
        Field('subtype', Reference())
    ]


class Record(ComplexNode):
    __type__ = 'Record'
    __name_field__ = 'type_name'
    __value_field__ = 'type_name'

    _fields = ComplexNode._fields + [
        Field('type_name'),
        Field('input_type'),
        Field('type_full_name'),
        Field('keys', TypedList(RecordKey)),
        Field('implements', TypedList(Reference)),
        Field('reducible_to_key')
    ]


class AbstractRecord(ComplexNode):
    __type__ = 'AbstractRecord'
    __name_field__ = 'name'
    __value_field__ = 'name'

    _fields = Record._fields + [
        Field('implementations', TypedList(Reference)),
        Field('default_descendant', Reference()),
        Field('full_name'),
        Field('name')
    ]


class Selection(ComplexNode):
    __type__ = 'Selection'
    __name_field__ = 'name'
    __value_field__ = 'name'

    _fields = ComplexNode._fields + [
        Field('values', TypedList(SelectionValue)),
        Field('name'),
        Field('full_name')
    ]

    def include_in_format(self):
        return self.name.find('TYPE') == -1


# all acceptable input_type
registered_nodes = {
    'Record': Record,
    'AbstractRecord': AbstractRecord,
    'Selection': Selection,
    'String': String,
    'Double': Double,
    'Integer': Integer,
    'FileName': FileName,
    'Bool': Bool,
    'Array': Array
}
