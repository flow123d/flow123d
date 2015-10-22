# encoding: utf-8
# author:   Jan Hybs
import cgi
from ist.nodes import Integer, String, DescriptionNode, ISTNode, ComplexNode
from ist.utils.htmltree import htmltree


class HTMLItemFormatter(htmltree):
    """
    Simple formatter class
    """
    def __init__(self, cls):
        super(HTMLItemFormatter, self).__init__('section', cls)


class HTMLUniversal(HTMLItemFormatter):
    """
    More abstract formatter class for strings, booleans, and other simple elements
    """
    def __init__(self):
        super(HTMLUniversal, self).__init__(cls='simple-element')

    def _start_format_as_child(self, self_object, record_key, record):
        self.h(record_key.key, record.type_name)

    def _format_as_child(self, self_object, record_key, record):
        raise Exception('Not implemented yet')

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
        method = self.format_rules.get(self_default.type, None)
        if method:
            return method(self_default, record_key, record)

        return HTMLRecordKeyDefault.textlangle_format(self_default, record_key, record)

    def textlangle_format(self, self_default, record_key, record):
        self.html.info('Default value: ')
        with self.html.open('span', attrib={ 'class': 'item-value chevron header-info skew' }):
            if len(str(self_default.value)):
                self.html.span(self_default.value)
            else:
                self.html.span(self_default.type)

        return self.html.current()

    def raw_format(self, self_default, record_key, record):
        self.html.info('Default value: ')
        with self.html.open('span', attrib={ 'class': 'item-value chevron skew' }):
            if len(str(self_default.value)):
                self.html.span(self_default.value)
            else:
                self.html.span(self_default.type)
        return self.html.current()


class HTMLInteger(HTMLUniversal):
    """
    Class representing int
    """
    def _format_as_child(self, self_int, record_key, record):
        with self.open('div', attrib={ 'class': 'item-key-value' }):
            with self.open('span', attrib={ 'class': 'item-key' }):
                self.info('Integer ')
            with self.open('span', attrib={ 'class': 'item-' }):
                self.span(str(self_int.range))


class HTMLDouble(HTMLUniversal):
    """
    Class representing double
    """
    def _format_as_child(self, self_double, record_key, record):
        with self.open('div', attrib={ 'class': 'item-key-value' }):
            with self.open('span', attrib={ 'class': 'item-key' }):
                self.info('Double ')
            with self.open('span', attrib={ 'class': 'item-value' }):
                self.span(str(self_double.range))


class HTMLBool(HTMLUniversal):
    """
    Class representing boolean
    """
    def _format_as_child(self, self_bool, record_key, record):
        with self.open('div', attrib={ 'class': 'item-key-value' }):
            with self.open('span', attrib={ 'class': 'item-key' }):
                self.span('Bool (generic)')

            self.tag('br')
            HTMLRecordKeyDefault(self).format_as_child(record_key.default, record_key, record)


class HTMLString(HTMLUniversal):
    """
    Class representing string
    """
    def _format_as_child(self, self_fn, record_key, record):
        with self.open('div', attrib={ 'class': 'item-key-value' }):
            with self.open('span', attrib={ 'class': 'item-key' }):
                self.span('String (generic)')


class HTMLFileName(HTMLUniversal):
    """
    Class representing filename type
    """
    def _format_as_child(self, self_fn, record_key, record):
        with self.open('div', attrib={ 'class': 'item-key-value' }):
            with self.open('span', attrib={ 'class': 'item-key' }):
                self.info('file name')
            with self.open('span', attrib={ 'class': 'item-value chevron' }):
                self.span(self_fn.file_mode)


class HTMLArray(HTMLUniversal):
    """
    Class representing Array structure
    """
    def _format_as_child(self, self_array, record_key, record):
        subtype = self_array.subtype.get_reference()
        with self.open('div', attrib={ 'class': 'item-key-value' }):
            with self.open('span', attrib={ 'class': 'item-key' }):
                self.info(' Array')
                if str(self_array.range):
                    self.info(' ')
                    self.span(cgi.escape(str(self_array.range)))
                self.info(' of ')
                self.span(subtype.get_type())

            if issubclass(subtype.__class__, ComplexNode):
                with self.open('span', attrib={ 'class': 'item-value chevron' }):
                    self.link(subtype.get_name())

            self.tag('br')
            HTMLRecordKeyDefault(self).format_as_child(record_key.default, record_key, record)

            #
            # if type(subtype) == Integer:
            # self.span('Array of {subtype} {subrange}'.format(
            # range=self_array.range, subtype=subtype.input_type,
            # subrange=subtype.range))
            # else:
            # self.span('Array{range} of {subtype}'.format(
            # range=' ' + str(self_array.range) if not self_array.range.is_pointless() else '',
            # subtype=subtype.input_type))
            #
            # if type(subtype) == String:
            #         self.span(' (generic)')
            #
            #     if issubclass(subtype.__class__, DescriptionNode):
            #         self.span(': ')
            #         self.link(subtype.get('type_name', 'name'))  # TODO href
            #     else:
            #         # no link
            #         pass


class HTMLSelection(HTMLItemFormatter):
    """
    Class representing Selection node in IST
    """
    def __init__(self):
        super(HTMLSelection, self).__init__(cls='main-section selection hidden')

    def format_as_child(self, self_selection, record_key, record):
        self.root.attrib['class'] = 'child-selection'
        with self.open('div', attrib={ 'class': 'item-key-value' }):
            with self.open('span', attrib={ 'class': 'item-key' }):
                self.info('Selection ')
            with self.open('span', attrib={ 'class': 'item-value chevron' }):
                if self_selection.include_in_format():
                    self.link(self_selection.get_name())
                else:
                    self.span(self_selection.get_name())

            self.tag('br')
            HTMLRecordKeyDefault(self).format_as_child(record_key.default, record_key, record)

        self.h(record_key.key, record.get_name())

        self.description(record_key.description)

    def format(self, selection):
        self.root.attrib['id'] = htmltree.chain_values(selection.get_name())
        self.root.attrib['data-name'] = htmltree.chain_values(selection.get_name())
        with self.open('header'):
            self.h2(selection.name)
            self.description(selection.description)

        with self.open('ul', attrib={ 'class': 'item-list' }):
            for selection_value in selection.values:
                with self.open('li'):
                    self.h3(selection_value.get_name())
                    self.description(selection_value.description)

        return self


class HTMLAbstractRecord(HTMLItemFormatter):
    """
    Class representing AbstractRecord node in IST
    """
    def __init__(self):
        super(HTMLAbstractRecord, self).__init__(cls='main-section abstract-record hidden')

    def format_as_child(self, abstract_record, record_key, record):
        self.root.attrib['class'] = 'child-abstract-record'
        with self.open('div', attrib={ 'class': 'item-key-value' }):
            with self.open('span', attrib={ 'class': 'item-key' }):
                self.info('abstract type ')
            with self.open('span', attrib={ 'class': 'item-value chevron' }):
                if abstract_record.include_in_format():
                    self.link(abstract_record.get_name())
                else:
                    self.span(abstract_record.get_name())

            self.tag('br')
            HTMLRecordKeyDefault(self).format_as_child(record_key.default, record_key, record)

        self.h(record_key.key, record.get_name())
        self.description(record_key.description)

    def format(self, abstract_record):
        self.root.attrib['id'] = htmltree.chain_values(abstract_record.get_name())
        self.root.attrib['data-name'] = htmltree.chain_values(abstract_record.get_name())
        with self.open('header'):
            self.h2(abstract_record.name)

            if abstract_record.default_descendant:
                reference = abstract_record.default_descendant.get_reference()
                self.italic('Default descendant ')
                self.link(reference.get_name())

            self.description(abstract_record.description)

        self.italic('Implemented by:')
        with self.open('ul', attrib={ 'class': 'item-list' }):
            for descendant in abstract_record.implementations:
                reference = descendant.get_reference()
                with self.open('li'):
                    self.link(reference.type_name)
                    self.info(' - ')
                    self.span(reference.description)


class HTMLRecord(HTMLItemFormatter):
    """
    Class representing record node in IST
    """
    def __init__(self):
        super(HTMLRecord, self).__init__(cls='main-section record hidden')

    def format_as_child(self, self_record, record_key, record):
        self.root.attrib['class'] = 'child-record'
        with self.open('div', attrib={ 'class': 'item-key-value' }):
            with self.open('span', attrib={ 'class': 'item-key' }):
                self.info('Record ')
            with self.open('span', attrib={ 'class': 'item-value chevron' }):
                if self_record.include_in_format():
                    self.link(self_record.get_name())
                else:
                    self.span(self_record.get_name())

            self.tag('br')
            HTMLRecordKeyDefault(self).format_as_child(record_key.default, record_key, record)

        self.h(record_key.key, record.get_name())

        self.description(record_key.description)

    def format(self, record):
        self.root.attrib['id'] = htmltree.chain_values(record.get_name())
        self.root.attrib['data-name'] = htmltree.chain_values(record.get_name())
        reference_list = record.implements
        with self.open('header'):
            self.h2(record.type_name)

            if reference_list:
                for reference in reference_list:
                    self.italic('implements abstract type: ')
                    self.link(reference.get_reference().name)

            if record.reducible_to_key:
                if reference_list:
                    self.tag('br')
                self.italic('constructible from key: ')
                self.link(record.reducible_to_key, ns=record.type_name)

            self.description(record.description)

        with self.open('ul', attrib={ 'class': 'item-list' }):
            for record_key in record.keys:
                if not record_key.include_in_format():
                    continue
                with self.open('li'):
                    fmt = HTMLFormatter.get_formatter_for(record_key)
                    fmt.format(record_key, record)
                    self.add(fmt.current())


class HTMLRecordKey(HTMLItemFormatter):
    """
    Class representing one record key
    """
    def __init__(self):
        super(HTMLRecordKey, self).__init__(cls='record-key')

    def format(self, record_key, record):
        reference = record_key.type.get_reference()

        # try to grab formatter and format type and default value based on reference type
        try:
            fmt = HTMLFormatter.get_formatter_for(reference)
            fmt.format_as_child(reference, record_key, record)
            self.add(fmt.current())
        except Exception as e:
            print ' <<Missing formatter for {}>>'.format(type(reference))
            print e


class HTMLFormatter(object):
    """
    Class which performs formatting
    """
    formatters = {
        'Record': HTMLRecord,
        'RecordKey': HTMLRecordKey,
        'AbstractRecord': HTMLAbstractRecord,
        'String': HTMLString,
        'Selection': HTMLSelection,
        'Array': HTMLArray,
        'Integer': HTMLInteger,
        'Double': HTMLDouble,
        'Bool': HTMLBool,
        'FileName': HTMLFileName,
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

        for item in items:
            # format only IST nodes
            if issubclass(item.__class__, ISTNode):

                # do no format certain objects
                if not item.include_in_format():
                    continue

                try:
                    fmt = HTMLFormatter.get_formatter_for(item)
                    fmt.format(item)
                    html.add(fmt.current())
                except Exception as e:
                    # print e
                    continue

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

        html.bold('Abstract record ')
        HTMLFormatter._add_items(sorted_items, html, 'AbstractRecord', reverse=False)

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
        prev_name = ''
        with html.open('ul', attrib={ 'class': 'nav-bar' }):
            for item in items:
                if issubclass(item.__class__, ComplexNode):
                    # do no format certain objects
                    if not item.include_in_format():
                        continue

                    if type and not item.input_type == type:
                        continue

                    if prev_name == item.get_name():
                        continue

                    prev_name = item.get_name()

                    with html.open('li', attrib={'data-name': item.get_name()}):
                        with html.open('a', '', html.generate_href(item.get_name())):
                            if reverse:
                                html.span(item.get_name())
                                html.span(item.get_type()[0], attrib={ 'class': 'shortcut-r' })
                            else:
                                html.span(item.get_type()[0], attrib={ 'class': 'shortcut' })
                                html.span(item.get_name())

    @staticmethod
    def __cmp(a, b):
        try:
            name_a = a.get_name()
        except:
            name_a = None

        try:
            name_b = b.get_name()
        except:
            name_b = None

        if name_a > name_b:
            return 1

        if name_a < name_b:
            return -1
        return 0
