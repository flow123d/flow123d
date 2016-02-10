#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

import cgi

import xml.etree.ElementTree as ET
import re
from ist.formatters.markdown2html import markdown2html


class htmltree(object):
    """
    Helper class for generating html documents
    """
    m2h = markdown2html()

    def __init__(self, tag_name='div', cls='', *args, **kwargs):
        self.attrib = { 'class': cls } if cls else { }
        self.tag_name = tag_name
        self.root = ET.Element(tag_name, self.attrib)
        self.counter = 0
        self.roots = [self.root]

    def tag(self, tag_name, value='', attrib={ }, no_escape=True, **kwargs):
        """
        Method will create and append element with tag based on tag_name value
        :param tag_name: tag name
        :param value: tag value
        :param attrib: optional attribute dict
        :param no_escape: whether to escape value
        :return: element
        """
        attrib_copy = { }
        if 'cls' in kwargs:
            attrib_copy['class'] = kwargs.pop('cls')

        attrib_copy.update(attrib)
        attrib_copy.update(kwargs)

        element = ET.Element(tag_name, attrib_copy)
        element.text = cgi.escape(value) if not no_escape else value
        self.current().append(element)
        return element

    def current(self):
        """
        :return: current top
        """
        return self.roots[self.counter]

    def add(self, element):
        """
        :return: appends element to top
        """
        return self.current().append(element)

    def item_list_title(self, item_field, add_link=False, add_id=True):
        """
        Method creates and appends header with level
        If subtitle is given href to title will consist of subtitle

        :param title: main header title
        :param subtitle: optional subtitle
        :param level: h level
        :param hide_subtitle: hide subtitle section
        :return:
        """
        if add_link:
            with self.open('a', attrib={'href': '#'+item_field.href_id}):
                self.item_list_title(item_field, add_link=False, add_id=add_id)
        else:
            attrib = {'id': item_field.href_id} if add_id else {}
            with self.open('h3', attrib=attrib):
                self.span(item_field.get('name', 'key', 'href_name'))
        return self

    def main_section_title(self, item, attrib={ }, **kwargs):
        """
        Method creates level 2 header also with "on the side" href with icon
          to this href
        :type item: ist.nodes.TypeSelection or ist.nodes.TypeRecord or ist.nodes.TypeAbstract
        """
        # with self.open('span', '', { 'class': 'pull-right side-anchor' }):
        #     href_attrib = { 'href': '#' + item.href_id }
        #     href_attrib.update({ 'title': 'Permalink to this section' })
        #     with self.open('a', '', href_attrib):
        #         self.span(' ', { 'class': 'glyphicon glyphicon-link', 'aria-hidden': 'true' })
        with self.open('h2'):
            self.tag('a', item.href_name, attrib={ 'href': '#' + item.href_id })
        # self.tag('h2', item.href_name, attrib, **kwargs)

    def mark_as_obsolete(self, element):
        """
        :type item: ist.nodes.TypeSelection or ist.nodes.TypeRecord or ist.nodes.TypeAbstract
        """
        with self.open('div', cls='obsolete'):
            self.tag('p', 'Obsolete', cls='obsolete-title')
            self.description(element.attributes.obsolete)

    def add_clear(self):
        self.div('', cls='clear')

    def h2(self, value='', attrib={ }, **kwargs):
        """
        Method creates level 3 header
        :param value: header title
        :param attrib: optional attribute
        :return: element
        """
        # attrib.update(self.generate_id(value))
        self.tag('h2', value, attrib, **kwargs)

    def h3(self, value='', attrib={ }, **kwargs):
        """
        Method creates level 3 header
        :param value: header title
        :param attrib: optional attribute
        :return: element
        """
        return self.tag('h3', value, attrib, **kwargs)

    def h4(self, value='', attrib={ }, **kwargs):
        """
        Method creates level 4 header
        :param value: header title
        :param attrib: optional attribute
        :return: element
        """
        return self.tag('h4', value, attrib, **kwargs)

    def h5(self, value='', attrib={ }, **kwargs):
        """
        Method creates level 5 header
        :param value: header title
        :param attrib: optional attribute
        :return: element
        """
        return self.tag('h5', value, attrib, **kwargs)

    def h6(self, value='', attrib={ }, **kwargs):
        """
        Method creates level 6 header
        :param value: header title
        :param attrib: optional attribute
        :return: element
        """
        return self.tag('h6', value, attrib, **kwargs)

    def ul(self, value='', attrib={ }, **kwargs):
        """
        Method creates ul element
        :param value: ul optional text
        :param attrib: optional attribute
        :return: element
        """
        return self.tag('ul', value, attrib, **kwargs)

    def ol(self, value='', attrib={ }, **kwargs):
        """
        Method creates ol element
        :param value: ol optional text
        :param attrib: optional attribute
        :return: element
        """
        return self.tag('ol', value, attrib, **kwargs)

    def span(self, value='', attrib={ }, **kwargs):
        """
        Method creates span element
        :param value: span optional text
        :param attrib: optional attribute
        :return: element
        """
        return self.tag('span', value, attrib, **kwargs)

    def info(self, value='', attrib={ "class": 'leading-text' }, **kwargs):
        """
        Method creates info element
        :param value: span optional text
        :param attrib: optional attribute, default has class of 'leading-text'
        :return: element
        """
        return self.tag('span', value, attrib, **kwargs)

    def div(self, value='', attrib={ }, **kwargs):
        """
        Method creates div element
        :param value: div optional text
        :param attrib: optional attribute
        :return: element
        """
        return self.tag('div', value, attrib, **kwargs)

    def bold(self, value='', attrib={ }, **kwargs):
        """
        Method creates bold element
        :param value: bold optional text
        :param attrib: optional attribute
        :return: element
        """
        return self.tag('strong', value, attrib, **kwargs)

    def italic(self, value='', attrib={ }, **kwargs):
        """
        Method creates em element
        :param value: em optional text
        :param attrib: optional attribute
        :return: element
        """
        return self.tag('em', value, attrib, **kwargs)

    def li(self, value='', attrib={ }, **kwargs):
        """
        Method creates li element
        :param value: li optional text
        :param attrib: optional attribute
        :return: element
        """
        return self.tag('li', value, attrib, **kwargs)

    def link(self, target, text='', ns=''):
        """
        Method creates link element based on given attributes
        :param target: machine link name
        :param text: link title
        :param ns: namespace for link
        """
        return self.tag('a', text if text else target, attrib=self.generate_href(target, ns))

    def link_to_main(self, item, text=None):
        """
        :type item: ist.nodes.TypeSelection or ist.nodes.TypeRecord or ist.nodes.TypeAbstract
        """
        # return self.tag('a', item.name + '(' + item.id + ')', self.generate_href(item.id))
        return self.tag('a', text or item.href_name, {'href': '#'+item.href_id})

    def open(self, tag_name, value='', attrib={ }, **kwargs):
        """
        Method opens current element, shifts current top.
          When adding new elements, current top is this newly created
        :param tag_name: tag name
        :param value: optional text
        :param attrib: optional attribs
        :return: self
        """
        element = self.tag(tag_name, value, attrib, **kwargs)
        self.roots.append(element)
        return self

    def description(self, value):
        """
        method will append description element to top
        element has predefined class
        :param value:
        :return:
        """
        # self.tag('span', 'DESCRIPTION', attrib={ 'class': 'desc-title' })
        if not value:
            return self.tag('div', 'no description provided', cls='description no-description')

        value = self.m2h.parse(value, reduce_to_tree=True)
        value.attrib['class'] = 'description'
        self.add(value)
        # return self.tag('div', value, { 'class': 'description' }, no_escape=True)

    def __enter__(self):
        """
        Enter the runtime context related to this object.
        :return:
        """
        self.counter += 1
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        """
        Exit the runtime context related to this object.
        :param exception_type:
        :param exception_value:
        :param traceback:
        :return:
        """
        # add debug info
        if exception_type:
            return False

        self.counter -= 1
        self.roots.pop()
        return self

    def dump(self):
        """
        Debug html dump
        :return:
        """
        return ET.tostring(self.root, method='html')

    def __repr__(self):
        return '<htmltree object>'

    def style(self, location):
        """
        Method add css link style
        :param location: css file location relative to server
        :return: element
        """
        self.tag('link', '', rel='stylesheet', type='text/css', media='screen', href=location)

    def script(self, location):
        """
        Method add css link script
        :param location: css file location relative to server
        :return: element
        """
        self.tag('script', '', attrib={ }, type='text/javascript', src=location)

    def id(self, id):
        """
        Method sets id to root element
        :param id: id value
        """
        self.root.attrib['id'] = id

    @staticmethod
    def secure(value):
        """
        Method for securing input value
        :param value: value to be secured
        :return: safe value consisting of numbers and alphabet chars and _ char
        """
        value = re.sub(r'\W+', '-', value)
        value = re.sub(r'-$', '', value)
        return value

    @staticmethod
    def chain_values(value, sub_value=''):
        """
        Method for chaining given values
        Used in href and id creation
        :param value: main value
        :param sub_value: optional sub value
        :return: secured chained value: "values-subvalue"  or "value" in lowercase
        """
        chain = value if not sub_value else sub_value + '-' + value
        return htmltree.secure(chain).lower()