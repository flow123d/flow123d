# encoding: utf-8
# author:   Jan Hybs
from __future__ import absolute_import

import cgi
from ist.globals import Globals
from ist.utils.htmltree import htmltree
from ist.base import InputType, NotImplementedException

from utils.logger import Logger


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


class HTMLSelection(HTMLItemFormatter):
    """
    Class representing Selection node in IST
    """

    def __init__(self):
        super(HTMLSelection, self).__init__(cls='main-section selection hidden')

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
            refs = selection.get_references()
            if refs:
                with self.open('div', cls='references'):
                    self.info('used in: ')
                    self.tag('br')
                    for ref in refs:
                        with self.open('span'):
                            self.link_to_main(ref)
                            if refs.index(ref) != len(refs)-1:
                                self.info(', ')
            self.main_section_title(selection)
            self.description(selection.description)

        if selection.values:
            self.italic('Values', attrib={'class': 'section-list'})
            with self.open('ul', attrib={'class': 'item-list'}):
                for selection_value in selection.values:
                    with self.open('li'):
                        self.item_list_title(selection_value)
                        self.description(selection_value.description)

        if selection.attributes.parameters:
            self.italic('Parameters', attrib={'class': 'section-list'})
            with self.open('ul', attrib={'class': 'item-list'}):
                for param in selection.attributes.parameters:
                    reference = param.reference.get_reference()
                    with self.open('li'):
                        with self.open('section', attrib={'class': 'record-param'}):
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
        super(HTMLRecord, self).__init__(cls='main-section record hidden')

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
            refs = record.get_references()
            if refs:
                with self.open('div', cls='references'):
                    self.info('used in: ')
                    self.tag('br')
                    for ref in refs:
                        with self.open('span'):
                            self.link_to_main(ref)
                            if refs.index(ref) != len(refs)-1:
                                self.info(', ')
            self.main_section_title(record)

            if record.attributes.generic_type:
                with self.open('div'):
                    self.italic('Generic type: ')
                    self.link_to_main(record.attributes.generic_type.get_reference())

            if record.reducible_to_key:
                with self.open('div'):
                    self.italic('Constructible from key: ')
                    self.link_to_main(record.reducible_to_key)

            if record.implements:
                with self.open('div'):
                    self.italic('Implements abstract type: ')
                    with self.open('ul'):
                        for reference in record.implements:
                            with self.open('li'):
                                self.link_to_main(reference.get_reference())

            self.italic('Description', attrib={'class': 'section-list'})
            self.description(record.description)

        if record.keys:
            self.italic('Keys', attrib={'class': 'section-list'})
            with self.open('ul', attrib={'class': 'item-list'}):
                for record_key in record.keys:
                    if not record_key.include_in_format():
                        continue
                    with self.open('li'):
                        with self.open('section', attrib={'class': 'record-key'}):
                            self.item_list_title(record_key, add_link=True)
                            self.format_key(record_key)
                            self.description(record_key.description)
                            self.add_clear()
                            # fmt = HTMLFormatter.get_formatter_for(record_key)
                            # fmt.format(record_key, record)
                            # self.add(fmt.current())

        if record.attributes.parameters:
            self.italic('Parameters', attrib={'class': 'section-list'})
            with self.open('ul', attrib={'class': 'item-list'}):
                for param in record.attributes.parameters:
                    reference = param.reference.get_reference()
                    with self.open('li'):
                        with self.open('section', attrib={'class': 'record-param'}):
                            self.h3(param.name)
                            self.span(str(reference.input_type))
                            self.info(' type of ')
                            self.link_to_main(reference)

    def format_key(self, record_key):
        """

        :type record_key: ist.extras.TypeRecordKey
        """
        if record_key.type.get_reference():
            reference = record_key.type.get_reference()
            input_type = reference.input_type

            if input_type == InputType.STRING:
                with self.open('ul', cls='side-info'):
                    with self.open('li'):
                        with self.open('div', attrib={'class': 'item-key-value'}):
                            with self.open('span', attrib={'class': 'item-key'}):
                                self.span('String (generic)')

            if input_type == InputType.BOOL:
                with self.open('ul', cls='side-info'):
                    with self.open('li'):
                        with self.open('div', attrib={'class': 'item-key-value'}):
                            with self.open('span', attrib={'class': 'item-key'}):
                                self.span('Bool (generic)')

            if input_type == InputType.DOUBLE:
                with self.open('ul', cls='side-info'):
                    with self.open('li'):
                        with self.open('div', attrib={'class': 'item-key-value'}):
                            with self.open('span', attrib={'class': 'item-key'}):
                                self.info('Double ')
                            with self.open('span', attrib={'class': 'item-value'}):
                                self.span(str(reference.range))

            if input_type == InputType.INTEGER:
                with self.open('ul', cls='side-info'):
                    with self.open('li'):
                        with self.open('div', attrib={'class': 'item-key-value'}):
                            with self.open('span', attrib={'class': 'item-key'}):
                                self.info('Integer ')
                            with self.open('span', attrib={'class': 'item-value'}):
                                self.span(str(reference.range))

            if input_type == InputType.FILENAME:
                with self.open('ul', cls='side-info'):
                    with self.open('li'):
                        with self.open('div', attrib={'class': 'item-key-value'}):
                            with self.open('span', attrib={'class': 'item-key'}):
                                self.info('file name')
                            with self.open('span', attrib={'class': 'item-value chevron'}):
                                self.span(reference.file_mode)

            if input_type == InputType.ARRAY:
                subtype = reference.subtype.get_reference()
                with self.open('ul', cls='side-info'):
                    with self.open('li'):
                        with self.open('span', attrib={'class': 'item-key'}):
                            self.info(' Array')
                            if str(reference.range):
                                self.info(' ')
                                self.span(cgi.escape(str(reference.range)))

                    with self.open('li'):
                        if not subtype.input_type == InputType.MAIN_TYPE:
                            self.info(' of ')
                            self.span(str(subtype.input_type))
                        else:
                            # self.tag('br')
                            if not subtype.attributes.generic_type:
                                self.info(' of ')
                                with self.open('span', attrib={'class': 'item-value chevron'}):
                                    self.link_to_main(subtype)
                        self.add_genericity(subtype)

            if input_type == InputType.MAIN_TYPE:
                with self.open('ul', cls='side-info'):
                    with self.open('li'):
                        if not reference.attributes.generic_type:
                            with self.open('span', attrib={'class': 'item-key'}):
                                self.span(str(reference.input_type))
                            with self.open('span', attrib={'class': 'item-value chevron'}):
                                self.link_to_main(reference)
                        else:
                            with self.open('span', attrib={'class': 'item-value chevron'}):
                                self.link_to_main(reference, text='instance')
                            self.info(' of ')
                            with self.open('span', attrib={'class': 'item-key'}):
                                self.span(str(reference.input_type))
                        self.add_genericity(reference)

    def add_genericity(self, item):
        if not item.attributes.generic_type:
            return
        generic_class = item.attributes.generic_type.get_reference()
        generic_impl = item
        with self.open('li'):
            self.info(' generic type')
            with self.open('span', attrib={'class': 'item-value chevron'}):
                self.link_to_main(generic_class)

        for param in generic_impl.attributes.parameters:
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


class HTMLAbstractRecord(HTMLItemFormatter):
    """
    Class representing AbstractRecord node in IST
    """

    def __init__(self):
        super(HTMLAbstractRecord, self).__init__(cls='main-section abstract-record hidden')

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
            refs = abstract_record.get_references()
            if refs:
                with self.open('div', cls='references'):
                    self.info('used in: ')
                    self.tag('br')
                    for ref in refs:
                        with self.open('span'):
                            self.link_to_main(ref)
                            if refs.index(ref) != len(refs)-1:
                                self.info(', ')

            self.main_section_title(abstract_record)

            if abstract_record.default_descendant:
                reference = abstract_record.default_descendant.get_reference()
                with self.open('div'):
                    self.italic('Default descendant: ')
                    self.link_to_main(reference)

            if abstract_record.attributes.generic_type:
                with self.open('div'):
                    self.italic('Generic type: ')
                    self.link_to_main(abstract_record.attributes.generic_type.get_reference())

            self.italic('Description', attrib={'class': 'section-list'})
            self.description(abstract_record.description)

        if abstract_record.implementations:
            self.italic('Implementations', attrib={'class': 'section-list'})
            with self.open('ul', attrib={'class': 'item-list'}):
                for descendant in abstract_record.implementations:
                    reference = descendant.get_reference()
                    with self.open('li'):
                        with self.open('section', attrib={'class': 'record-param'}):
                            with self.open('h3'):
                                self.link_to_main(reference)
                            self.span(reference.description)

        if abstract_record.attributes.parameters:
            self.italic('Parameters', attrib={'class': 'section-list'})
            with self.open('ul', attrib={'class': 'item-list'}):
                for param in abstract_record.attributes.parameters:
                    reference = param.reference.get_reference()
                    with self.open('li'):
                        with self.open('section', attrib={'class': 'record-param'}):
                            self.h3(param.name)
                            self.span(str(reference.input_type))
                            self.info(' type of ')
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
        html.id('input-reference')
        Logger.instance().info('Processing items...')

        for item in items:

            # do no format certain objects
            if not item.include_in_format():
                Logger.instance().info(' - item skipped: %s' % str(item))
                continue

            Logger.instance().info(' - formatting item: %s' % str(item))
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
        html = htmltree('div')

        sorted_items = sorted(items, cmp=HTMLFormatter.__cmp)

        html.bold('Records ')
        HTMLFormatter._add_items(sorted_items, html, 'Record', reverse=False)

        html.bold('Abstract records ')
        HTMLFormatter._add_items(sorted_items, html, 'Abstract', reverse=False)

        html.bold('Selections ')
        HTMLFormatter._add_items(sorted_items, html, 'Selection', reverse=False)

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
    def _add_items(items, html, type=None, reverse=False):
        """
        :type items: list[ist.nodes.TypeSelection]
        """
        prev_name = ''
        with html.open('ul', attrib={'class': 'nav-bar'}):
            for item in items:
                if item.input_type == InputType.MAIN_TYPE:
                    # do no format certain objects
                    if not item.include_in_format():
                        continue

                    if type and not item.input_type == type:
                        continue

                    if prev_name == item.name:
                        continue

                    prev_name = item.name

                    with html.open('li', attrib={'data-name': item.name}):
                        with html.open('a', '', {'href': '#' + item.href_id, 'class': 'item_' + item.href_id}):
                            if reverse:
                                html.span(item.href_name)
                                html.span(str(item.input_type)[0], attrib={'class': 'shortcut-r'})
                            else:
                                html.span(str(item.input_type)[0], attrib={'class': 'shortcut'})
                                html.span(item.href_name)

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
