#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

from __future__ import absolute_import
from ist.base import InputType
from utils.logger import Logger
from ist.nodes import TypeInteger, TypeString
from ist.utils.texlist import texlist


class LatexItemFormatter(object):
    def __init__(self, tag_name=None):
        self.tag_name = tag_name

    def format(self, element):
        raise NotImplementedError('Method format not implemented {}'.format(self.__class__.__name__))

    def format_as_child(self, *args):
        raise NotImplementedError('Method format_as_child not implemented {}'.format(self.__class__.__name__))


class LatexSelection(LatexItemFormatter):
    def __init__(self):
        super(LatexSelection, self).__init__('SelectionType')

    def format_as_child(self, self_selection, record_key, record):
        """
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        :type self_selection: ist.nodes.TypeSelection
        """
        tex = texlist()
        tex.KeyItem()
        with tex:
            tex.hyperB(record_key.key, record.name)

        with tex:
            tex.append('selection: ')
            if self_selection.include_in_format():
                tex.Alink(self_selection.name)
            else:
                tex.append(tex.escape(self_selection.name))

        tex.add_s(record_key.default.value)
        tex.add()
        tex.add_description_field(record_key.description)

        return tex

    def format(self, selection):
        """
        :type selection: ist.nodes.TypeSelection
        """
        tex = texlist(self.tag_name)

        with tex.element():
            with tex:
                tex.hyperB(selection.name)
            tex.add_description_field(selection.description)

            for selection_value in selection.values:
                tex.newline()
                tex.KeyItem(selection_value.name, selection_value.description)
            tex.newline()

        return tex


class LatexRecordKey(LatexItemFormatter):
    def __init__(self):
        super(LatexRecordKey, self).__init__('KeyItem')

    def format(self, record_key, record):
        """
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        """
        tex = texlist()
        reference = record_key.type.get_reference()

        # try to grab formatter and format type and default value based on reference type
        fmt = LatexFormatter.get_formatter_for(reference)
        if not fmt:
            Logger.instance().info(' <<Missing formatter for {}>>'.format(type(reference)))
        else:
            tex.extend(fmt.format_as_child(reference, record_key, record))

        return tex


class LatexAbstractRecord(LatexItemFormatter):
    def __init__(self):
        super(LatexAbstractRecord, self).__init__('AbstractType')

    def format_as_child(self, abstract_record, record_key, record):
        """
        :type abstract_record: ist.nodes.TypeAbstract
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        """
        tex = texlist()
        tex.KeyItem()
        with tex:
            tex.hyperB(record_key.key, record.name)

        with tex:
            tex.append('abstract type: ')
            tex.Alink(abstract_record.name)

        tex.extend(
            LatexRecordKeyDefault().format_as_child(record_key.default, record_key, record)
        )
        tex.add()
        tex.add_description_field(record_key.description)

        return tex

    def format(self, abstract_record):
        """
        :type abstract_record: ist.nodes.TypeAbstract
        """
        tex = texlist(self.tag_name)
        with tex.element():
            with tex:
                tex.hyperB(abstract_record.name)

            if abstract_record.default_descendant:
                reference = abstract_record.default_descendant.get_reference()
                with tex:
                    tex.Alink(reference.name)
                with tex:
                    tex.AddDoc(abstract_record.name)
            else:
                tex.add()
                tex.add()
            tex.add_description_field(abstract_record.description)

            for descendant in abstract_record.implementations:
                tex.newline()
                tex.tag('Descendant')
                with tex:
                    tex.Alink(descendant.get_reference().name)
            tex.newline()

        return tex


class LatexString(LatexItemFormatter):
    def __init__(self):
        super(LatexString, self).__init__('String')

    def format_as_child(self, self_string, record_key, record):
        """
        :type self_fn: ist.nodes.TypeString
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        """
        tex = texlist()
        tex.KeyItem()
        with tex:
            tex.hyperB(record_key.key, record.name)
        tex.add('String (generic)')

        with tex:
            tex.textlangle(record_key.default.type)
        tex.add()
        tex.add_description_field(record_key.description)

        return tex


class LatexRecord(LatexItemFormatter):
    def __init__(self):
        super(LatexRecord, self).__init__('RecordType')

    def format_as_child(self, self_record, record_key, record):
        """
        :type self_record: ist.nodes.TypeRecord
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        """
        tex = texlist()
        tex.KeyItem()
        with tex:
            tex.hyperB(record_key.key, record.name)

        with tex:
            tex.append('record: ')
            tex.Alink(self_record.name)

        tex.extend(
            LatexRecordKeyDefault().format_as_child(record_key.default, record_key, record)
        )
        tex.add()
        tex.add_description_field(record_key.description)

        return tex

    def format(self, record):
        """
        :type record: ist.nodes.TypeRecord
        """
        tex = texlist(self.tag_name)
        reference_list = record.implements

        with tex.element():
            with tex:
                tex.hyperB(record.name)

            # TODO what if multiple inheritance? list
            if reference_list:
                with tex:
                    for reference in reference_list:
                        tex.Alink(reference.get_reference().name)
            else:
                tex.add()

            if record.reducible_to_key:
                with tex:
                    tex.Alink(record.reducible_to_key, record.name)
            else:
                tex.add()

            # hyperlink into hand written text
            # LATER it can removed since is not used anymore
            tex.add()
            tex.add_description_field(record.description)

            # record keys
            for record_key in record.keys:
                tex.newline()
                fmt = LatexFormatter.get_formatter_for(record_key)
                tex.extend(fmt.format(record_key, record))
            tex.newline()

        return tex


class LatexRecordKeyDefault(LatexItemFormatter):
    def __init__(self):
        super(LatexRecordKeyDefault, self).__init__('')

        self.format_rules = {
            'value at read time': self.raw_format,
            'value at declaration': self.textlangle_format,
            'optional': self.textlangle_format,
            'obligatory': self.textlangle_format
        }

    def format_as_child(self, self_default, record_key, record):
        """
        :type self_default: ist.extras.TypeRecordKeyDefault
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        """
        method = self.format_rules.get(self_default.type, None)
        if method:
            return method(self_default, record_key, record)

        return LatexRecordKeyDefault.textlangle_format(self_default, record_key, record)

    def textlangle_format(self, self_default, record_key, record):
        """
        :type self: object
        :type self_default: ist.extras.TypeRecordKeyDefault
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        """
        tex = texlist()
        with tex:
            tex.textlangle(self_default.value)
        return tex

    def raw_format(self, self_default, record_key, record):
        """
        :type self_default: ist.extras.TypeRecordKeyDefault
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        """
        tex = texlist()
        tex.add_s(self_default.value)
        return tex


class LatexUniversal(LatexItemFormatter):
    def __init__(self):
        super(LatexUniversal, self).__init__('')

    def _start_format_as_child(self, self_object, record_key, record):
        """
        :type self_record: object
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        """
        tex = texlist()
        tex.KeyItem()
        with tex:
            tex.hyperB(record_key.key, record.name)
        return tex

    def _end_format_as_child(self, self_object, record_key, record):
        """
        :type self_record: object
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        """
        tex = texlist()

        tex.extend(
            LatexRecordKeyDefault().format_as_child(record_key.default, record_key, record)
        )
        tex.add()
        tex.add_description_field(record_key.description)

        return tex

    def _format_as_child(self, self_object, record_key, record):
        """
        :type self_record: object
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        """
        raise Exception('Not implemented yet')

    def format_as_child(self, self_object, record_key, record):
        """
        :type self_record: object
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        """
        tex = texlist()
        tex.extend(self._start_format_as_child(self_object, record_key, record))
        tex.extend(self._format_as_child(self_object, record_key, record))
        tex.extend(self._end_format_as_child(self_object, record_key, record))
        return tex


class LatexArray(LatexUniversal):
    def _format_as_child(self, self_array, record_key, record):
        """
        :type self_array: ist.nodes.TypeArray
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        """
        subtype = self_array.subtype.get_reference()
        tex = texlist()
        with tex:
            if type(subtype) == TypeInteger:
                tex.append('Array of {subtype} {subrange}'.format(
                    range=self_array.range, subtype=subtype.input_type,
                    subrange=subtype.range))
            else:
                tex.append('Array{range} of {subtype}'.format(
                    range=' ' + str(self_array.range) if not self_array.range.is_pointless() else '',
                    subtype=subtype.input_type))

            if type(subtype) == TypeString:
                tex.append(' (generic)')

            if subtype.input_type == InputType.MAIN_TYPE:
                tex.append(': ')
                tex.Alink(subtype.name)
            else:
                # no link
                pass

        return tex


class LatexInteger(LatexUniversal):
    def _format_as_child(self, self_int, record_key, record):
        """
        :type self_int: ist.nodes.TypeInteger
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        """
        tex = texlist()
        tex.add('Integer' + str(self_int.range))
        return tex

    def _end_format_as_child(self, self_object, record_key, record):
        """
        :type self_int: ist.nodes.TypeInteger
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        """
        tex = texlist()
        tex.extend(
            LatexRecordKeyDefault().format_as_child(record_key.default, record_key, record)
        )
        tex.add()
        tex.add_description_field(record_key.description)
        return tex


class LatexDouble(LatexUniversal):
    def _format_as_child(self, self_double, record_key, record):
        """
        :type self_double: ist.nodes.TypeDouble
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        """
        tex = texlist()
        tex.add('Double' + str(self_double.range))
        return tex

    def _end_format_as_child(self, self_object, record_key, record):
        """
        :type self_double: ist.nodes.TypeDouble
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        """
        tex = texlist()
        tex.extend(
            LatexRecordKeyDefault().format_as_child(record_key.default, record_key, record)
        )
        tex.add()
        tex.add_description_field(record_key.description)
        return tex


class LatexBool(LatexUniversal):
    def _format_as_child(self, self_bool, record_key, record):
        """
        :type self_bool: ist.nodes.TypeBool
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        """
        tex = texlist()
        tex.add('Bool')
        return tex

    def _end_format_as_child(self, self_object, record_key, record):
        """
        :type self_bool: ist.nodes.TypeBool
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        """
        tex = texlist()
        tex.add(record_key.default.value)  # todo LatexRecordKeyDefault
        tex.add()
        tex.add_description_field(record_key.description)

        return tex


class LatexFileName(LatexUniversal):
    def _format_as_child(self, self_fn, record_key, record):
        """
        :type self_fn: ist.nodes.TypeFilename
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        """
        tex = texlist()
        tex.add(self_fn.file_mode + ' file name')
        return tex


class LatexFormatter(object):
    formatters = {
        'TypeRecord': LatexRecord,
        'TypeRecordKey': LatexRecordKey,
        'TypeAbstract': LatexAbstractRecord,
        'TypeString': LatexString,
        'TypeSelection': LatexSelection,
        'TypeArray': LatexArray,
        'TypeInteger': LatexInteger,
        'TypeDouble': LatexDouble,
        'TypeBool': LatexBool,
        'TypeFilename': LatexFileName
    }

    @staticmethod
    def format(items):
        tex = texlist()

        Logger.instance().info('Processing items...')
        for item in items:
            # do no format certain objects
            if not item.include_in_format():
                Logger.instance().info(' - item skipped: %s' % str(item))
                continue

            Logger.instance().info(' - formatting item: %s' % str(item))
            fmt = LatexFormatter.get_formatter_for(item)
            if fmt is not None:
                tex.extend(fmt.format(item))
                tex.newline()
                tex.newline()

        return tex

    @staticmethod
    def get_formatter_for(o):
        cls = LatexFormatter.formatters.get(o.__class__.__name__, None)
        if cls is None:
            return None
        return cls()