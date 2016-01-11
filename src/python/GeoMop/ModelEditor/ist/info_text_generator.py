"""
Module generates html documentation from IST.
"""
import json
from .utils.htmltree import htmltree
from copy import copy

__author__ = 'Tomas Krizek'


class InfoTextGenerator:
    """
    Generates info_text for `DataNode`.
    """
    _input_types = {}

    @classmethod
    def init(cls, json_text):
        """Initializes the class with format information."""
        data = json.loads(json_text, encoding="utf-8")
        for item in data:
            if 'id' in item:
                cls._input_types[item['id']] = item

    @classmethod
    def get_info_text(cls, record_id=None, selected_key=None, abstract_id=None, selected_item=None,
                      context=None):
        """Generates HTML documentation for `record_id` with `selected_key`.

        The first key is selected by default. If `abstract_id` is specified, a documentation for
        a parent abstract record is generated. If the selected key type is of type `Selection`,
        then `selected_item` will be selected.

        `context` is a dictionary containing query data for home, back and forward buttons.
        """
        cls.record_id = record_id
        cls.selected_key = selected_key
        cls.abstract_id = abstract_id
        cls.selected_item = selected_item

        html = htmltree('html')
        html_body = htmltree('body')

        with html_body.open('div', cls='container-fluid fill'):
            if abstract_id in cls._input_types:
                html_body.add(cls._generate_abstract_record(abstract_id))
            if record_id in cls._input_types:
                html_body.add(cls._generate_record(record_id, selected_key, selected_item))
            if abstract_id not in cls._input_types and record_id not in cls._input_types:
                with html_body.open('section', cls='row record'):
                    html_body.description('')

        if context is not None and len(context.get('home', {})) > 0:
            # are we on different page from home?
            home = context['home']
            check_fields = ['record_id', 'abstract_id', 'selected_key', 'selected_item']
            displaying_home = True
            for field in check_fields:
                if locals().get(field) != home.get(field):
                    displaying_home = False
                    break

            if not displaying_home:  # show navigation bar
                with html_body.open('div', cls='navigation-panel'):
                    if len(context['back']) > 0:
                        args = copy(context['back'][-1])
                        args.update({'direction': 'back'})
                        href = cls.generate_href(**args)
                        class_ = 'back'
                    else:
                        href = '#'
                        class_ = 'back not-active'
                    html_body.tag('a', '', attrib={'href': href, 'class': class_})

                    if len(context['forward']) > 0:
                        args = copy(context['forward'][-1])
                        args.update({'direction': 'forward'})
                        href = cls.generate_href(**args)
                        class_ = 'forward'
                    else:
                        href = '#'
                        class_ = 'forward not-active'
                    html_body.tag('a', '', attrib={'href': href, 'class': class_})

                    if 'home' in context:
                        args = copy(context['home'])
                        args.update({'direction': 'home'})
                        href = cls.generate_href(**args)
                        class_ = 'home'
                    html_body.tag('a', '', attrib={'href': href, 'class': class_})

        html_head = htmltree('head')
        html_head.style('bootstrap.min.css')
        html_head.style('katex.min.css')
        html_head.style('main.css')
        html_head.script('jquery-2.1.3.min.js')
        html_head.script('bootstrap.min.js')

        html_body.script('katex.min.js')
        html_body.script('main.js')

        html.add(html_head.current())
        html.add(html_body.current())
        return r'<!DOCTYPE html>' + html.dump()

    @classmethod
    def _generate_abstract_record(cls, abstract_id):
        """Generates documentation for top-level abstract record."""
        section = htmltree('section', cls='row abstract-record')
        type_ = cls._input_types[abstract_id]

        with section.open('header'):
            section.tag('h2', type_.get('name', ''))
            section.description(type_.get('description', ''))

        cls._generate_implementation_list(section, type_)

        return section.current()

    @classmethod
    def _generate_record(cls, record_id, selected_key=None, selected_item=None):
        """Generates documentation for record."""
        section = htmltree('section', cls='row record')
        type_ = cls._input_types[record_id]
        selected_key_type = None

        with section.open('header'):
            section.tag('h2', type_.get('type_name', ''))
            section.description(type_.get('description', ''))

        with section.open('div', cls='key-list col-md-4 col-sm-4 col-xs-4'):
            with section.open('div', cls='item-list'):
                if 'keys' in type_:
                    keys = [key for key in type_['keys'] if key['key'] not in ['TYPE']]
                    if selected_key is None and len(keys) > 0:
                        selected_key = keys[0]['key']
                    for key in keys:
                        if key['key'] == selected_key:
                            selected_key_type = key
                            cls_ = 'selected'
                        else:
                            cls_ = ''
                        href = cls.generate_href(
                            record_id=record_id,
                            selected_key=key['key'],
                            abstract_id=cls.abstract_id,
                            selected_item=cls.selected_item
                        )
                        # change color/font for different key types
                        if 'default' in key:
                            cls_ += ' key-type-' + key['default']['type'].replace(' ', '-')
                        section.tag('a', key['key'], attrib={'class': cls_, 'href': href})

        if selected_key_type is not None:
            with section.open('div', cls='key-description col-md-4 col-sm-4 col-xs-4'):
                with section.open('header'):
                    section.tag('h3', selected_key)
                    if 'default' in selected_key_type and 'value' in selected_key_type['default']:
                        with section.open('div', cls='small'):
                            section.info('Default value: ')
                            section.tag('span', selected_key_type['default']['value'],
                                        cls='chevron skew')
                    section.description(selected_key_type.get('description', ''))

            key_type = cls._input_types.get(selected_key_type['type'])

            array_div = None
            while key_type is not None and key_type.get('input_type') == 'Array':
                if array_div is None:
                    array_div = htmltree('div', cls='small')
                array_div.info('Array ')
                range_ = NumberRange(key_type)
                array_div.tag('span', str(range_))
                array_div.info(' of ')
                array_div.tag('br')
                key_type = cls._input_types.get(key_type['subtype'])

            if key_type is not None and 'input_type' in key_type:
                cls_ = 'key-type col-md-4 col-sm-4 col-xs-4 '
                if key_type['input_type'] == 'Record':
                    cls_ += 'record'
                    section.add(cls._generate_key_type_record(key_type, cls_=cls_,
                                                              prepend=array_div))
                elif key_type['input_type'] == 'AbstractRecord':
                    cls_ += 'abstract-record'
                    section.add(cls._generate_key_type_abstract_record(key_type, cls_=cls_,
                                                                       prepend=array_div))
                else:
                    cls_ += 'scalar'
                    if key_type['input_type'] == 'Selection':
                        section.add(cls._generate_key_type_selection(key_type, selected_item,
                                                                     cls_=cls_, prepend=array_div))
                    else:
                        section.add(cls._generate_key_type_scalar(key_type, cls_=cls_,
                                                                  prepend=array_div))

        return section.current()

    @classmethod
    def _generate_key_type_scalar(cls, key_type, cls_='', prepend=None):
        """Generates documentation for scalar key type."""
        div = htmltree('div', cls=cls_)
        if prepend is not None:
            div.add(prepend.current())
        with div.open('header'):
            div.tag('h2', key_type.get('name', ''))
            if key_type.get('input_type') in ['Integer', 'Double']:
                range_ = NumberRange(key_type)
                div.tag('span', str(range_), cls='small')
            elif key_type.get('input_type') == 'FileName':
                with div.open('div', cls='small'):
                    div.info('File mode: ')
                    div.span(key_type.get('file_mode', 'unknown'), cls='chevron skew')
            else:
                div.description('')
        return div.current()

    @classmethod
    def _generate_key_type_selection(cls, key_type, selected_item=None, cls_='', prepend=None):
        """Generates documentation for scalar key type."""
        div = htmltree('div', cls=cls_)
        if prepend is not None:
            div.add(prepend.current())
        with div.open('header'):
            div.tag('h2', key_type.get('name', ''))
            div.description(key_type.get('description', ''))

        values = {value['name']: value['description'] for value in key_type.get('values', [])}
        if selected_item in values:
            div.description('<span class="leading-text">{0}: </span>{1}'
                            .format(selected_item, values[selected_item]))

        with div.open('div', cls='item-list'):
            for name in values:
                if name == selected_item:
                    cls_ = 'selected'
                else:
                    cls_ = ''
                href = cls.generate_href(
                    record_id=cls.record_id,
                    selected_key=cls.selected_key,
                    abstract_id=cls.abstract_id,
                    selected_item=name
                )
                div.tag('a', name, attrib={'class': cls_, 'href': href})

        return div.current()

    @classmethod
    def _generate_key_type_record(cls, key_type, cls_='', prepend=None):
        """Generates documentation for record key type."""
        div = htmltree('div', cls=cls_)
        if prepend is not None:
            div.add(prepend.current())
        with div.open('header'):
            href = cls.generate_href(
                record_id=key_type.get('id')
            )
            div.tag('a', key_type.get('type_name', ''), attrib={'class': 'h2', 'href': href})
            if 'reducible_to_key' in key_type:
                with div.open('div', cls='small'):
                    div.info('Constructible from key: ')
                    div.span(key_type['reducible_to_key'], cls='chevron skew')
        div.description(key_type.get('description', ''))
        return div.current()

    @classmethod
    def _generate_key_type_abstract_record(cls, key_type, cls_='', prepend=None):
        """Generates documentation for abstract record key type."""
        div = htmltree('div', cls=cls_)
        if prepend is not None:
            div.add(prepend.current())
        with div.open('header'):
            href = cls.generate_href(
                abstract_id=key_type.get('id')
            )
            div.tag('a', key_type.get('name', ''), attrib={'class': 'h2', 'href': href})
            div.description(key_type.get('description', ''))
        cls._generate_implementation_list(div, key_type)
        return div.current()

    @classmethod
    def _generate_implementation_list(cls, tree, type_):
        """Generates implementation list for abstract record."""
        with tree.open('ul', cls='implementation-list'):
            for implementation_id in type_.get('implementations', []):
                implementation_type = cls._input_types.get(implementation_id)
                if implementation_type is None:
                    continue
                with tree.open('li'):
                    name = implementation_type.get('type_name')
                    if not name:
                        name = implementation_type.get('name', '')
                    href = cls.generate_href(
                        record_id=implementation_id,
                        abstract_id=type_['id']
                    )
                    tree.tag('a', name, attrib={'href': href})
                    tree.info(' - ')
                    tree.span(implementation_type.get('description', ''))

    @staticmethod
    def generate_href(record_id=None, selected_key=None, abstract_id=None,
                      selected_item=None, direction=None):
        """Generates href link from data."""
        # pylint: disable=unused-argument
        parts = []
        for name in ['record_id', 'selected_key', 'abstract_id', 'selected_item']:
            value = locals().get(name)
            if value is not None:
                parts.append("{0}={1}".format(name, value))

        if direction is not None:
            parts.append("{0}=1".format(direction))

        if not parts:
            return '#'
        return '?' + '&'.join(parts)


class NumberRange:
    """
    Class representing simple number range
    """

    def __init__(self, input_type):
        self.min = self.max = ''
        if 'range' in input_type:
            self.min = input_type['range'][0]
            self.max = input_type['range'][1]

    replacements = {
        '2147483647': 'INT32 MAX',
        '4294967295': 'UINT32 MAX',
        '-2147483647': 'INT32 MIN',
        '1.79769e+308': '+inf',
        '-1.79769e+308': '-inf',
        '': 'unknown range'
    }

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

    def __repr__(self):
        """Returns string representation."""
        return self._format()

