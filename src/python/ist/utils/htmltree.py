# encoding: utf-8
# author:   Jan Hybs
import cgi

import xml.etree.ElementTree as ET
import re
from ist.formatters.markdown2html import markdown2html


class htmltree(object):
    m2h = markdown2html()

    def __init__(self, tag_name='div', cls='', *args, **kwargs):
        self.attrib = { 'class': cls } if cls else { }
        self.tag_name = tag_name
        self.root = ET.Element(tag_name, self.attrib)
        self.counter = 0
        self.roots = [self.root]

    def tag(self, tag_name, value='', attrib={ }, no_escape=True):
        element = ET.Element(tag_name, attrib)
        element.text = cgi.escape(value) if not no_escape else value
        self.current().append(element)
        return element

    def current(self):
        return self.roots[self.counter]

    def add(self, element):
        return self.current().append(element)

    def h(self, title, subtitle='', level='h3', hide_sub=True):
        if subtitle:
            with self.open(level, '', self.generate_id(title, subtitle)):
                if not hide_sub:
                    with self.open('small'):
                        self.tag('a', subtitle + '', self.generate_href(subtitle))
                        self.span('::')
                self.span(title)
            return self
        with self.open(level, ''):
            with self.open('a', '', self.generate_href(title)):
                self.span(title)
        return self

    def h1(self, value='', attrib={ }):
        return self.tag('h1', value, attrib)

    def h2(self, value='', attrib={ }):
        with self.open('span', '', { 'class': 'pull-right side-anchor' }):
            href_attrib = self.generate_href(value)
            href_attrib.update({ 'title': 'Permalink to this section' })
            with self.open('a', '', href_attrib):
                self.span(' ', { 'class': 'glyphicon glyphicon-link', 'aria-hidden': 'true' })

        # attrib.update(self.generate_id(value))
        self.tag('h2', value, attrib)

    def h3(self, value='', attrib={ }):
        return self.tag('h3', value, attrib)

    def h4(self, value='', attrib={ }):
        return self.tag('h4', value, attrib)

    def h5(self, value='', attrib={ }):
        return self.tag('h5', value, attrib)

    def h6(self, value='', attrib={ }):
        return self.tag('h6', value, attrib)

    def ul(self, value='', attrib={ }):
        return self.tag('ul', value, attrib)

    def ol(self, value='', attrib={ }):
        return self.tag('ol', value, attrib)

    def span(self, value='', attrib={ }):
        return self.tag('span', value, attrib)

    def info(self, value='', attrib={"class": 'leading-text'}):
        return self.tag('span', value, attrib)

    def div(self, value='', attrib={ }):
        return self.tag('div', value, attrib)

    def bold(self, value='', attrib={ }):
        return self.tag('strong', value, attrib)

    def italic(self, value='', attrib={ }):
        return self.tag('em', value, attrib)

    def li(self, value='', attrib={ }):
        return self.tag('li', value, attrib)

    def link(self, target, text='', ns=''):
        return self.tag('a', text if text else target, self.generate_href(target, ns))

    def open(self, tag_name, value='', attrib={ }):
        element = self.tag(tag_name, value, attrib)
        self.roots.append(element)
        return self

    def description(self, value):
        if not value:
            return self.tag('div', 'no description provided', { 'class': 'description no-description' })

        value = self.m2h.parse(value, reduce_to_tree=True)
        value.attrib['class'] = 'description'
        self.current().append(value)
        # return self.tag('div', value, { 'class': 'description' }, no_escape=True)


    def __enter__(self):
        self.counter += 1
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.counter -= 1
        self.roots.pop()
        return self

    def dump(self):
        return ET.tostring(self.root, method='html')

    def __repr__(self):
        return '<htmltree object>'

    def style(self, location):
        self.tag('link', '', { 'rel': 'stylesheet', 'type': 'text/css', 'media': 'screen', 'href': location })

    def script(self, location):
        self.tag('script', '', { 'type': 'text/javascript', 'src': location })

    def id(self, id):
        self.root.attrib['id'] = id

    def generate_id(self, value, sub_value=''):
        return { 'id': htmltree.chain_values(value, sub_value) }

    def generate_href(self, value, sub_value=''):
        return { 'href': '#' + htmltree.chain_values(value, sub_value) }

    @staticmethod
    def secure(value):
        value = re.sub(r'\W+', '-', value)
        value = re.sub(r'-$', '', value)
        return value

    @staticmethod
    def chain_values(value, sub_value=''):
        chain = value if not sub_value else sub_value + '-' + value
        return htmltree.secure(chain).lower()