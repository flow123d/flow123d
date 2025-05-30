#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs


import cgi
import html
from py123d.ist.utils.htmltree import htmltree
from py123d.ist.base import InputType, NotImplementedException
from py123d.utils.logger import Logger


class HTMLItemFormatter(htmltree):
    """
    Simple formatter class
    """

    def __init__(self, cls):
        super(HTMLItemFormatter, self).__init__('section', cls)

    def format_as_child(self, *args, **kwargs):
        raise NotImplementedException('Not implemented yet')

    def format(self, *args, **kwargs):
        raise NotImplementedException('Not implemented yet')

    def add_header_right_side(self, parsable):
        with self.openc('div', 'info'):
            refs = list(set(parsable.references))
            if refs:
                self.info('used in: ')
                with self.openc('div'):
                    for ref in refs:
                        with self.open('span'):
                            self.link_to_main(ref)
                            if refs.index(ref) != len(refs) - 1:
                                self.info(', ')
            fers = list(set(parsable.secnerefer))
            if fers:
                self.info('referenced by:')
                with self.openc('div'):
                    for fer in fers:
                        with self.open('span'):
                            self.link_to_main(fer)
                            if fers.index(fer) != len(fers) - 1:
                                self.info(', ')


class HTMLSelection(HTMLItemFormatter):
    """
    Class representing Selection node in IST
    """

    def __init__(self):
        super(HTMLSelection, self).__init__(cls='s')

    def format_as_child(self, self_selection, record_key, record):
        """
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        :type self_selection: ist.nodes.TypeSelection
        """
        self.root.attrib['class'] = 'child-selection'
        with self.open('div', attrib={'class': 'item-key-value'}):
            with self.open('span', attrib={'class': 'item-key'}):
                self.info('Selection ')
            with self.open('span', attrib={'class': 'item-value chevron'}):
                if self_selection.include_in_format():
                    self.link_to_main(self_selection)
                else:
                    self.span(self_selection.name)

            self.tag('br')
            HTMLRecordKeyDefault(self).format_as_child(record_key.default, record_key, record)

        self.item_list_title(record_key)
        self.description(record_key.description)

    def format(self, selection):
        """
        :type selection: ist.nodes.TypeSelection
        """
        self.root.attrib['id'] = selection.href_id
        self.root.attrib['data-name'] = htmltree.secure(selection.name)

        if selection.attributes.obsolete:
            self.root.attrib['data-obsolete'] = '1'
            self.mark_as_obsolete(selection)

        with self.open('header'):
            self.add_header_right_side(selection)
            self.main_section_title(selection)

            with self.openc('div', 'details'):
                with self.openc('div', 'section'):
                    self.info('Description')
                    self.description(selection.description)

        if selection.values:
            self.spanc('section-list', 'Values')
            with self.open('ul', attrib={'class': 'item-list'}):
                for selection_value in selection.values:
                    with self.open('li'):
                        with self.openc('div', 'item-header'):
                            self.item_list_title(selection_value, add_link=True)
                        self.description(selection_value.description)

        if selection.parameters:
            self.spanc('section-list', 'Parameters')
            with self.open('ul', attrib={'class': 'item-list'}):
                for param in selection.parameters:
                    reference = param.reference.get_reference()
                    with self.open('li'):
                        self.h3(param.name)
                        self.span(str(reference.input_type))
                        self.info(' type of ')
                        self.link_to_main(reference)

        return self


class HTMLRecord(HTMLItemFormatter):
    """
    Class representing record node in IST
    """

    def __init__(self):
        super(HTMLRecord, self).__init__(cls='r')

    def format_as_child(self, self_record, record_key, record):
        """
        :type self_record: ist.nodes.TypeRecord
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        """
        self.root.attrib['class'] = 'child-record'
        with self.open('div', attrib={'class': 'item-key-value'}):
            with self.open('span', attrib={'class': 'item-key'}):
                self.info('Record ')
            with self.open('span', attrib={'class': 'item-value chevron'}):
                if self_record.include_in_format():
                    self.link_to_main(self_record)
                else:
                    self.span(self_record.name)

            self.tag('br')
            HTMLRecordKeyDefault(self).format_as_child(record_key.default, record_key, record)

        self.item_list_title(record_key)

        self.description(record_key.description)

    def format(self, record):
        """
        :type record: ist.nodes.TypeRecord
        """
        self.root.attrib['id'] = record.href_id
        self.root.attrib['data-name'] = htmltree.secure(record.name)

        if record.attributes.obsolete:
            self.root.attrib['data-obsolete'] = '1'
            self.mark_as_obsolete(record)

        with self.open('header'):
            self.add_header_right_side(record)
            self.main_section_title(record)
            with self.openc('div', 'details'):
                if record.generic_type:
                    with self.openc('div', 'section'):
                        self.info('Generic type: ')
                        self.link_to_main(record.generic_type.get_reference())

                if record.reducible_to_key:
                    with self.openc('div', 'section'):
                        self.info('Constructible from key: ')
                        self.link_to_main(record.reducible_to_key)

                if record.implements:
                    with self.openc('div', 'section'):
                        self.info('Implements abstract type: ')
                        with self.open('ul'):
                            for reference in record.implements:
                                with self.open('li'):
                                    self.link_to_main(reference.get_reference())

                with self.openc('div', 'section'):
                    self.info('Description')
                    self.description(record.description)

        if record.keys:
            self.spanc('section-list', 'Keys')
            with self.openc('ul', 'item-list'):
                for record_key in record.keys:
                    if not record_key.include_in_format():
                        continue
                    with self.open('li'):
                        with self.open('div', attrib={'class': 'r-key'}):
                            link = self.format_key(record_key)
                            with self.openc('div', 'item-header'):
                                if link:
                                    self.item_list_title(link, add_link=True, add_id=False, text=record_key.href_name)
                                else:
                                    self.item_list_title(record_key, add_link=False)
                                self.format_key_default(record_key)

                            self.description(record_key.description)
                            self.add_clear()
                            # fmt = HTMLFormatter.get_formatter_for(record_key)
                            # fmt.format(record_key, record)
                            # self.add(fmt.current())

        if record.attributes.generic_parameters:
            self.spanc('section-list', 'Generic parameters')
            with self.open('ul', attrib={'class': 'item-list'}):
                for param_name in record.attributes.generic_parameters:
                    with self.open('li'):
                        with self.openc('div', 'item-header'):
                            self.h3(param_name)

        if record.parameters:
            self.spanc('section-list', 'Parameters')
            with self.open('ul', attrib={'class': 'item-list'}):
                for param in record.parameters:
                    reference = param.reference.get_reference()
                    with self.open('li'):
                        with self.openc('div', 'item-header'):
                            self.h3(param.name)
                            self.span(str(reference.input_type))
                            self.info(' type of ')
                            self.link_to_main(reference)

    def format_key_default(self, record_key):
        """
        :type record_key: ist.extras.TypeRecordKey
        """
        default = record_key.default
        if default.type in ('obligatory', 'optional'):
            self.spanc('chevron robot', str(default.type))
        else:
            self.spanc('chevron info', str(default.type))
            self.info(' = ')
            self.spanc('robot', str(default.value))

    def format_key(self, record_key):
        """
        :type record_key: ist.extras.TypeRecordKey
        """
        if record_key.type.get_reference():
            reference = record_key.type.get_reference()
            input_type = reference.input_type

            if input_type == InputType.STRING:
                with self.openc('ul', 'info'):
                    with self.open('li'):
                        with self.open('div', attrib={'class': 'item-key-value'}):
                            with self.open('span', attrib={'class': 'item-key'}):
                                self.span('String (generic)')

            if input_type == InputType.BOOL:
                with self.openc('ul', 'info'):
                    with self.open('li'):
                        with self.open('div', attrib={'class': 'item-key-value'}):
                            with self.open('span', attrib={'class': 'item-key'}):
                                self.span('Bool (generic)')

            if input_type == InputType.DOUBLE:
                with self.openc('ul', 'info'):
                    with self.open('li'):
                        with self.open('div', attrib={'class': 'item-key-value'}):
                            with self.open('span', attrib={'class': 'item-key'}):
                                self.info('Double ')
                            with self.open('span', attrib={'class': 'item-value'}):
                                self.span(str(reference.range))

            if input_type == InputType.INTEGER:
                with self.openc('ul', 'info'):
                    with self.open('li'):
                        with self.open('div', attrib={'class': 'item-key-value'}):
                            with self.open('span', attrib={'class': 'item-key'}):
                                self.info('Integer ')
                            with self.open('span', attrib={'class': 'item-value'}):
                                self.span(str(reference.range))

            if input_type == InputType.FILENAME:
                with self.openc('ul', 'info'):
                    with self.open('li'):
                        with self.open('div', attrib={'class': 'item-key-value'}):
                            with self.open('span', attrib={'class': 'item-key'}):
                                self.info('file name')
                            with self.open('span', attrib={'class': 'item-value chevron'}):
                                self.span(reference.file_mode)

            if input_type == InputType.ARRAY:
                subtype = reference.subtype.get_reference()
                with self.openc('ul', 'info'):
                    with self.open('li'):
                        with self.open('span', attrib={'class': 'item-key'}):
                            self.info(' Array')
                            if str(reference.range):
                                self.info(' ')
                                self.span(html.escape(str(reference.range)))

                    with self.open('li'):
                        if not subtype.input_type == InputType.MAIN_TYPE:
                            self.info(' of ')
                            self.span(str(subtype.input_type))
                        else:
                            if not subtype.has_generic_link():
                                self.info(' of ')
                                with self.open('span', attrib={'class': 'item-value chevron'}):
                                    self.link_to_main(subtype)
                                    return subtype
                            else:
                                return self.add_genericity(subtype)

            if input_type == InputType.MAIN_TYPE:
                with self.openc('ul', 'info'):
                    with self.open('li'):
                        if not reference.has_generic_link():
                            with self.open('span', attrib={'class': 'item-key'}):
                                self.span(str(reference.input_type))
                            with self.open('span', attrib={'class': 'item-value chevron'}):
                                self.link_to_main(reference)
                                return reference
                        else:
                            self.spanc('item-value', 'instance ')
                            self.span('of ')
                            # with self.open('span', attrib={'class': 'item-value chevron'}):
                            #     self.link_to_main(reference, text='instance')
                            # self.info(' of ')
                            with self.open('span', attrib={'class': 'item-key'}):
                                self.span(str(reference.input_type))
                            return self.add_genericity(reference)

    def add_genericity(self, item):
        if not item.has_generic_link():
            return
        generic_class = item.generic_type.get_reference()
        generic_impl = item
        with self.open('li'):
            self.info(' generic type')
            with self.open('span', attrib={'class': 'item-value chevron'}):
                self.link_to_main(generic_class)

        for param in generic_impl.parameters:
            # self.tag('br')
            with self.open('li'):
                self.info('parameter ')
                self.span(param.name)
                self.info(' = ')
                # with self.open('span', attrib={'class': 'item-value chevron'}):
                reference = param.reference.get_reference()
                if not reference.input_type == InputType.MAIN_TYPE:
                    self.span(str(reference.input_type))
                else:
                    self.link_to_main(reference)

        return generic_class


class HTMLTuple(HTMLRecord):

    def __init__(self):
        super(HTMLRecord, self).__init__(cls='t')


class HTMLAbstractRecord(HTMLItemFormatter):
    """
    Class representing AbstractRecord node in IST
    """

    def __init__(self):
        super(HTMLAbstractRecord, self).__init__(cls='a')

    def format_as_child(self, abstract_record, record_key, record):
        """
        :type abstract_record: ist.nodes.TypeAbstract
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        """
        self.root.attrib['class'] = 'child-abstract-record'
        with self.open('div', attrib={'class': 'item-key-value'}):
            with self.open('span', attrib={'class': 'item-key'}):
                self.info('abstract type ')
            with self.open('span', attrib={'class': 'item-value chevron'}):
                if abstract_record.include_in_format():
                    self.link_to_main(abstract_record)
                else:
                    self.span(abstract_record.name)

            self.tag('br')
            HTMLRecordKeyDefault(self).format_as_child(record_key.default, record_key, record)

        self.item_list_title(record_key)
        self.description(record_key.description)

    def format(self, abstract_record):
        """
        :type abstract_record: ist.nodes.TypeAbstract
        """
        self.root.attrib['id'] = abstract_record.href_id
        self.root.attrib['data-name'] = htmltree.secure(abstract_record.name)

        if abstract_record.attributes.obsolete:
            self.root.attrib['data-obsolete'] = '1'
            self.mark_as_obsolete(abstract_record)

        with self.open('header'):
            self.add_header_right_side(abstract_record)
            self.main_section_title(abstract_record)

            with self.openc('div', 'details'):

                if abstract_record.attributes.root_of_generic_subtree:
                    with self.openc('div', 'section'):
                        self.italic('This abstract is root of the generic tree')

                if abstract_record.default_descendant:
                    reference = abstract_record.default_descendant.get_reference()
                    with self.openc('div', 'section'):
                        self.info('Default descendant: ')
                        self.link_to_main(reference)

                if abstract_record.generic_type:
                    with self.openc('div', 'section'):
                        self.info('Generic type: ')
                        self.link_to_main(abstract_record.generic_type.get_reference())

                with self.openc('div', 'section'):
                    self.info('Description')
                    self.description(abstract_record.description)

        if abstract_record.implementations:
            self.spanc('section-list', 'Implementations')
            with self.open('ul', attrib={'class': 'item-list'}):
                for descendant in abstract_record.implementations:
                    reference = descendant.get_reference()
                    with self.open('li'):
                        with self.openc('div', 'item-header'):
                            with self.open('h3'):
                                self.link_to_main(reference)
                        self.description(reference.description)

        if abstract_record.attributes.generic_parameters:
            self.spanc('section-list', 'Generic parameters')
            with self.open('ul', attrib={'class': 'item-list'}):
                for param_name in abstract_record.attributes.generic_parameters:
                    with self.open('li'):
                        with self.openc('div', 'item-header'):
                            self.h3(param_name)

        if abstract_record.parameters:
            self.spanc('section-list', 'Parameters')
            with self.open('ul', attrib={'class': 'item-list'}):
                for param in abstract_record.parameters:
                    reference = param.reference.get_reference()
                    with self.open('li'):
                        with self.openc('div', 'item-header'):
                            self.h3(param.name)
                            self.span(' ' + str(reference.input_type))
                            self.info(' type of ')
                            with self.openc('span', 'chevron'):
                                self.link_to_main(reference)


class HTMLUniversal(HTMLItemFormatter):
    """
    More abstract formatter class for strings, booleans, and other simple elements
    """

    def __init__(self):
        super(HTMLUniversal, self).__init__(cls='simple-element')

    def _start_format_as_child(self, self_object, record_key, record):
        self.item_list_title(record_key)

    def _format_as_child(self, self_object, record_key, record):
        raise NotImplementedException('Not implemented yet')

    def _end_format_as_child(self, self_object, record_key, record):
        # HTMLRecordKeyDefault(self).format_as_child(record_key.default, record_key, record)
        self.description(record_key.description)

    def format_as_child(self, self_object, record_key, record):
        """
        Format this child in listing
        :param self_object:
        :param record_key:
        :param record:
        :return:
        """
        self._format_as_child(self_object, record_key, record)
        self._start_format_as_child(self_object, record_key, record)
        self._end_format_as_child(self_object, record_key, record)


class HTMLRecordKeyDefault(object):
    """
    Class representing default value in record key
    """

    def __init__(self, html):
        self.html = html if html is not None else htmltree('div', 'record-key-default')
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

        return HTMLRecordKeyDefault.textlangle_format(self_default, record_key, record)

    def textlangle_format(self, self_default, record_key, record):
        """
        :type self_default: ist.extras.TypeRecordKeyDefault
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey

        """
        self.html.info('Default value: ')
        with self.html.open('span', attrib={'class': 'item-value chevron header-info skew'}):
            if len(str(self_default.value)):
                self.html.span(str(self_default.value))
            else:
                self.html.span(str(self_default.type))

        return self.html.current()

    def raw_format(self, self_default, record_key, record):
        """
        :type self_default: ist.extras.TypeRecordKeyDefault
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey

        """
        self.html.info('Default value: ')
        with self.html.open('span', attrib={'class': 'item-value chevron skew'}):
            if len(str(self_default.value)):
                self.html.span(str(self_default.value))
            else:
                self.html.span(str(self_default.type))
        return self.html.current()


class HTMLInteger(HTMLUniversal):
    """
    Class representing int
    """

    def _format_as_child(self, self_int, record_key, record):
        """
        :type self_int: ist.nodes.TypeInteger
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        """
        with self.open('div', attrib={'class': 'item-key-value'}):
            with self.open('span', attrib={'class': 'item-key'}):
                self.info('Integer ')
            with self.open('span', attrib={'class': 'item-'}):
                self.span(str(self_int.range))


class HTMLDouble(HTMLUniversal):
    """
    Class representing double
    """

    def _format_as_child(self, self_double, record_key, record):
        """
        :type self_double: ist.nodes.TypeDouble
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        """
        with self.open('div', attrib={'class': 'item-key-value'}):
            with self.open('span', attrib={'class': 'item-key'}):
                self.info('Double ')
            with self.open('span', attrib={'class': 'item-value'}):
                self.span(str(self_double.range))


class HTMLBool(HTMLUniversal):
    """
    Class representing boolean
    """

    def _format_as_child(self, self_bool, record_key, record):
        """
        :type self_bool: ist.nodes.TypeBool
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        """
        with self.open('div', attrib={'class': 'item-key-value'}):
            with self.open('span', attrib={'class': 'item-key'}):
                self.span('Bool (generic)')

            self.tag('br')
            HTMLRecordKeyDefault(self).format_as_child(record_key.default, record_key, record)


class HTMLString(HTMLUniversal):
    """
    Class representing string
    """

    def _format_as_child(self, self_fn, record_key, record):
        """
        :type self_fn: ist.nodes.TypeString
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        """
        with self.open('div', attrib={'class': 'item-key-value'}):
            with self.open('span', attrib={'class': 'item-key'}):
                self.span('String (generic)')


class HTMLFileName(HTMLUniversal):
    """
    Class representing filename type
    """

    def _format_as_child(self, self_fn, record_key, record):
        """
        :type self_fn: ist.nodes.TypeFilename
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        """
        with self.open('div', attrib={'class': 'item-key-value'}):
            with self.open('span', attrib={'class': 'item-key'}):
                self.info('file name')
            with self.open('span', attrib={'class': 'item-value chevron'}):
                self.span(self_fn.file_mode)


class HTMLArray(HTMLUniversal):
    """
    Class representing Array structure
    """

    def _format_as_child(self, self_array, record_key, record):
        """
        :type self_array: ist.nodes.TypeArray
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        """
        subtype = self_array.subtype.get_reference()
        with self.open('div', attrib={'class': 'item-key-value'}):
            with self.open('span', attrib={'class': 'item-key'}):
                self.info(' Array')
                if str(self_array.range):
                    self.info(' ')
                    self.span(cgi.escape(str(self_array.range)))
                self.info(' of ')
                self.span(str(subtype.input_type))

            if subtype.input_type == InputType.MAIN_TYPE:
                with self.open('span', attrib={'class': 'item-value chevron'}):
                    self.link_to_main(subtype)

            self.tag('br')
            HTMLRecordKeyDefault(self).format_as_child(record_key.default, record_key, record)


class HTMLRecordKey(HTMLItemFormatter):
    """
    Class representing one record key
    """

    def __init__(self):
        super(HTMLRecordKey, self).__init__(cls='record-key')

    def format(self, record_key, record):
        """
        :type record: ist.nodes.TypeRecord
        :type record_key: ist.extras.TypeRecordKey
        """
        reference = record_key.type.get_reference()

        # try to grab formatter and format type and default value based on reference type
        try:
            fmt = HTMLFormatter.get_formatter_for(reference)
            fmt.format_as_child(reference, record_key, record)
            self.add(fmt.current())
        except NotImplementedException as e:
            Logger.instance().info(' <<Missing formatter for {}>>'.format(type(reference)))
            # raise e


class HTMLFormatter(object):
    """
    Class which performs formatting
    """
    formatters = {
        'TypeRecord': HTMLRecord,
        'TypeTuple': HTMLTuple,
        'TypeRecordKey': HTMLRecordKey,
        'TypeAbstract': HTMLAbstractRecord,
        'TypeString': HTMLString,
        'TypeSelection': HTMLSelection,
        'TypeArray': HTMLArray,
        'TypeInteger': HTMLInteger,
        'TypeDouble': HTMLDouble,
        'TypeBool': HTMLBool,
        'TypeFilename': HTMLFileName,
        '': HTMLUniversal
    }

    @staticmethod
    def format(items):
        """
        Formats given items to HTML format
        :param items: json items
        :return: html div element containing all given elements
        """
        html = htmltree('div')
        html.id('ist')
        Logger.instance().info('Processing items...')
        i = 0

        for item in items:
            i += 1
            if i > 1050:
                break
            # do no format certain objects
            if not item.include_in_format():
                Logger.instance().info('[SKIP] %s skipped\n' % str(item))
                continue

            Logger.instance().info('[ OK ] formatting item %s\n' % str(item))
            fmt = HTMLFormatter.get_formatter_for(item)

            fmt.format(item)
            html.add(fmt.current())

        return html

    @staticmethod
    def get_formatter_for(o):
        """
        Return formatter for given object
        :param o:
        :return: formatter
        """
        cls = HTMLFormatter.formatters.get(o.__class__.__name__, None)
        if cls is None:
            cls = HTMLFormatter.formatters.get('')
        return cls()

    @staticmethod
    def abc_navigation_bar(items):
        """
        Returns alphabet ordered links to given items
        :param items: all items
        :return: html div object
        """
        html = htmltree('div', 'menu-items')
        import functools
        sorted_items = sorted(items, key=functools.cmp_to_key(HTMLFormatter.__cmp))

        html.h3('Records')
        HTMLFormatter._add_items(sorted_items, html, 'Record', reverse=False, cls='r')

        html.h3('Tuples')
        HTMLFormatter._add_items(sorted_items, html, 'Tuple', reverse=False, cls='t')

        html.h3('Abstract records')
        HTMLFormatter._add_items(sorted_items, html, 'Abstract', reverse=False, cls='a')

        html.h3('Selections')
        HTMLFormatter._add_items(sorted_items, html, 'Selection', reverse=False, cls='s')

        return html

    @staticmethod
    def tree_navigation_bar(items):
        """
        Returns default ordered links to given items
        :param items: all items
        :return: html div object
        """
        html = htmltree('div')

        # sorted_items = sorted(items, cmp=HTMLFormatter.__cmp)
        html.bold('All items ')
        HTMLFormatter._add_items(items, html)

        return html

    @staticmethod
    def _add_items(items, html, type=None, reverse=False, cls=''):
        """
        :type items: list[ist.nodes.TypeSelection]
        """
        prev_name = ''
        with html.open('ul', attrib={'class': 'nav-bar ' + cls}):
            for item in items:
                if item.input_type == InputType.MAIN_TYPE:
                    # do no format certain objects
                    if not item.include_in_format():
                        continue

                    if type and not item.input_type == type:
                        continue

                    if prev_name == item.name and 0:
                        continue

                    prev_name = item.name

                    with html.open('li', attrib={'data-name': item.name}):
                        html.tag('a', item.href_name, {'href': '#' + item.href_id, 'class': 'item_' + item.href_id})

    @staticmethod
    def __cmp(a, b):
        try:
            name_a = a.name
        except:
            name_a = None

        try:
            name_b = b.name
        except:
            name_b = None

        if name_a > name_b:
            return 1

        if name_a < name_b:
            return -1
        return 0
