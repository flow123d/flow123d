# encoding: utf-8
# author:   Jan Hybs
from ist.globals import Globals
from ist.utils.utils import Field, TypedList, AttributeDict


class ISTNode(object):
    """
    Parent of all nodes in IST
    """
    # type of the object (similar to __class__.__name__)
    __type__ = ''

    # simple fields
    _fields = [
        Field('id'),
        Field('type_name'),
        Field('input_type'),
    ]

    def __init__(self):
        pass


    def parse(self, o={ }):
        """
        method parses given object
        :param o: object to parse
        :return: self
        """
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
            self.__setattr__(field.save_as, value)
            if field.name == 'id' and value is not None:
                Globals.items[value] = self

        return self

    def include_in_format(self):
        """
        :return: whether this ISTNode will appear in output specification
        """
        return True

    def get(self, *args):
        """
        Getter which will return given field value
        If not found other field passed as arguments are tested/returned
        If no field exists raise Exception
        :param args: field names
        :return: field value or raise
        """
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
        return '<{self.__class__.__name__}[{self.id}] {self._fields}>'.format(self=self)


class Reference(ISTNode):
    """
    Class representing reference object
    """
    __type__ = 'Reference'

    def __init__(self):
        self.ref_id = None


    def parse(self, o={ }):
        self.ref_id = o
        return self

    def get_reference(self):
        """
        :return: reference if exists otherwise None
        """
        return Globals.items.get(str(self.ref_id), None)

    def copy(self):
        return Reference()

    def __repr__(self):
        return "<Reference ({self.ref_id}) -> {ref}>".format(self=self, ref=self.get_reference())

    def __nonzero__(self):
        return self.get_reference() is not None


class NumberRange(ISTNode):
    """
    Class representing simple number range
    """
    __type__ = 'Range'

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
        Whether is information within this instance beneficial
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
        """
        return self._format() if self.always_visible or not self.is_pointless() else ''


class DoubleRange(NumberRange):
    """
    Class representing Double number range
    """
    __type__ = 'Range'

    def __init__(self, always_visible=False):
        super(DoubleRange, self).__init__(always_visible)

    def is_pointless(self):
        return self._format() == '[, ]'

    def copy(self):
        return DoubleRange(self.always_visible)


class Attribute(ISTNode):
    """
    Class representing attribute node
    """
    __type__ = 'Attribute'


class AttributeNode(ISTNode):
    """
    Class representing the most complex nodes (records selections abstract)
with description and supporting attributes
    """
    _fields = ISTNode._fields + [
        Field('attributes', AttributeDict(Attribute))
    ]


class SelectionValue(AttributeNode):
    """
    Class representing selection node
    """
    __type__ = 'Selection value'

    _fields = AttributeNode._fields + [
        Field('description'),
        Field('name')
    ]


class RecordKeyDefault(ISTNode):
    """
    Class representing default value in record key
    """
    __type__ = 'Defaults'

    _fields = [
        Field('type'),
        Field('value')
    ]


class RecordKey(ISTNode):
    """
    Class representing one record key
    """
    __type__ = 'Record key'
    __name_field__ = 'key'
    __value_field__ = 'default'

    _fields = [
        Field('key'),
        Field('type', Reference()),
        Field('default', RecordKeyDefault()),
        Field('description')
    ]

    def include_in_format(self):
        return self.key.find('TYPE') == -1


class String(AttributeNode):
    """
    Class representing string
    """
    __type__ = 'String'


class Double(AttributeNode):
    """
    Class representing double
    """
    __type__ = 'Double'
    __name_field__ = ''
    __value_field__ = 'range'

    _fields = AttributeNode._fields + [
        Field('range', DoubleRange())
    ]


class Integer(AttributeNode):
    """
    Class representing int
    """
    __type__ = 'Integer'
    __name_field__ = ''
    __value_field__ = 'range'

    _fields = AttributeNode._fields + [
        Field('range', NumberRange())
    ]


class FileName(AttributeNode):
    """
    Class representing filename type
    """
    __type__ = 'FileName'
    __name_field__ = ''
    __value_field__ = 'file_mode'

    _fields = AttributeNode._fields + [
        Field('file_mode')
    ]


class Bool(AttributeNode):
    """
    Class representing boolean
    """
    __type__ = 'Bool'

    _fields = AttributeNode._fields + [
    ]


class Array(AttributeNode):
    """
    Class representing Array structure
    """
    __type__ = 'Array'
    __name_field__ = ''
    __value_field__ = 'range'

    _fields = AttributeNode._fields + [
        Field('range', NumberRange(False)),
        Field('subtype', Reference())
    ]


class ComplexNode(AttributeNode):
    """
    Class for Record, Selection and AbstracRecord
    """
    __type__ = 'Complex type'

    _fields = AttributeNode._fields + [
        Field('description')
    ]

    def __repr__(self):
        return '<{self.__class__.__name__}[{self.id}] {self.type_name}>'.format(self=self)


class Record(ComplexNode):
    """
    Class representing record node in IST
    """
    __type__ = 'Record'

    _fields = ComplexNode._fields + [
        Field('keys', TypedList(RecordKey)),
        Field('implements', TypedList(Reference)),
        Field('reducible_to_key')
    ]


class AbstractRecord(ComplexNode):
    """
    Class representing AbstractRecord node in IST
    """
    __type__ = 'AbstractRecord'

    _fields = ComplexNode._fields + [
        Field('implementations', TypedList(Reference)),
        Field('default_descendant', Reference()),
    ]


class Selection(ComplexNode):
    """
    Class representing Selection node in IST
    """
    __type__ = 'Selection'

    _fields = ComplexNode._fields + [
        Field('values', TypedList(SelectionValue)),
    ]

    def include_in_format(self):
        return self.type_name.find('TYPE') == -1


class Parameter(AttributeNode):
    """
    Class representing Parameter node in IST
    """
    __type__ = 'Selection'

    _fields = AttributeNode._fields + [
    ]

# all acceptable input_type values
registered_nodes = {
    'Record': Record,
    'AbstractRecord': AbstractRecord,
    'Abstract': AbstractRecord,
    'Selection': Selection,
    'String': String,
    'Double': Double,
    'Integer': Integer,
    'FileName': FileName,
    'Bool': Bool,
    'Array': Array,
    'Parameter': Parameter
}
