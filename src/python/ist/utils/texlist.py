#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

from __future__ import absolute_import
import re
from ist.base import Unicode
from ist.formatters.html2latex import Html2Latex
from ist.formatters.markdown2html import markdown2html

class texlist(list):
    """
    Helper class for creating latex document
    """
    m2h = markdown2html()

    def __init__(self, name=''):
        super(list, self).__init__()
        self.name = name
        self.counter = 1

    def tag(self, field_name, *values):
        """
        Method adds \name{value1}{value2}
        :param field_name:
        :param values:
        :return: self
        """
        self.slash(field_name)
        for value in values:
            self.add_s(value)

        return self

    def KeyItem(self, name='', description=''):
        """
        Method add KeyItem tag
        :param name:
        :param description:
        :return:
        """
        self.slash('KeyItem')
        if name:
            self.add_s(name)
        if description:
            with self:
                self.add_description_field(description)

        return self

    def add(self, value=''):
        """
        Method add value surrounded with braces
        :param value:
        :return:
        """
        self.append("{" + str(value) + "}")
        return self

    def add_s(self, value=''):
        '''
        Add field with value escaped
        '''
        self.append("{" + self.name_mode(value) + "}")
        return self

    def add_d(self, value=''):
        '''
        Add field with value underscores replaced by dashed
        '''
        self.append("{" + self.name_mode(value) + "}")
        return self

    def open(self):
        """
        Add open brace
        """
        self.append('{')
        return self

    def close(self):
        """
        Add close brace
        """
        self.append('}')
        return self

    def hyperB(self, value, ns='IT::'):
        """
        Add HyperB element
        :param value: id
        :param ns: optional namespace
        :return: self
        """
        if __debug__:
            self.tag('hyperB', self.name_mode((ns if ns.endswith('::') else ns + '::') + value))
            self.add(self.plain_mode((ns if ns.endswith('::') else ns + '::') + value))
        else:
            self.tag('hyperB', self.name_mode((ns if ns.endswith('::') else ns + '::') + value))
            self.add(self.plain_mode(value))
        return self

    def slash(self, value=''):
        """
        Add \value
        :param value:
        :return: self
        """
        self.append('\\')
        if value:
            self.append(value)

        return self

    def Alink(self, url, ns="IT::", text=None):
        """
        Method adds Alink section
        :param url: Alink url
        :param ns: optional namespace
        :param text: optional text (otherwise url is used)
        :return: self
        """
        ns = str(ns)
        url = str(url)
        ns = ns if ns.endswith('::') else ns + '::'
        ns = ns if url.find('::') == -1 else ''

        if __debug__:
            self.tag('Alink', self.name_mode(ns + url))
            self.add(self.plain_mode(ns + (url if text is None else text)))
        else:
            self.tag('Alink', self.name_mode(ns + url))
            self.add(self.plain_mode(url if text is None else text))

        return self

    def AddDoc(self, value):
        """
        Adds \AddDoc{value}
        :return: self
        """
        self.slash('AddDoc')
        self.add(self.plain_mode(value))

    def textlangle(self, value, namespace='\\it '):
        """
        Add < value > with optional italic
        :param value:
        :param namespace:
        :return:
        """
        self.slash('textlangle')
        self.add(self.plain_mode(namespace + str(value) + ' ').lower())
        self.slash('textrangle')

        return self

    def newline(self):
        """
        Adds newline
        :return: self
        """
        self.append('\n')
        return self

    def element(self):
        """
        Resets counter
        :return:
        """
        self.counter = 0
        return self

    def open_element(self, name):
        """
        Opens current element by name
        \begin{name}
        :param name: element name
        :return: self
        """
        self.tag('begin', name)
        return self

    def close_element(self, name):
        """
        Closes current element by name
        \end{name}
        :param name: element name
        :return: self
        """
        self.tag('end', name)
        return self

    def __enter__(self):
        """
        Enter the runtime context related to this object.
        :return:
        """
        if self.counter == 0:
            self.open_element(self.name)
        else:
            self.open()
        self.counter += 1
        return self

    def __exit__(self, exception_type, exception_value, tb):
        """
        Exit the runtime context related to this object.
        :param exception_type:
        :param exception_value:
        :param tb:
        :return:
        """
        # add debug info
        if exception_type:
            return False

        self.counter -= 1
        if self.counter == 0:
            self.close_element(self.name)
        else:
            self.close()
        return self

    def add_description_field(self, value):
        """
        Adds complex description field
        :param value: string value with markdown support
        :return: self
        """
        self.add(self.description(value))

    def description(self, value):
        """
        Creates complex description field
        :param value: string value with markdown support
        :return: self
        """
        # return self.escape (value.strip ().replace ('\n', '\\\\'))
        html = self.m2h.parse(str(value), True)
        latex = Html2Latex(html)
        result = latex.to_latex()
        desc_result = list()
        for r in result:
            if r.startswith('{$') and r.endswith('$}'):
                desc_result.append(self.equation_mode(r))
            else:
                desc_result.append(self.plain_mode(r))
        return ''.join(desc_result)


    @staticmethod
    def equation_mode(value):
        """
        Method will ensure that value will be valid equation in latex
        :type value: str or list
        :param value: value tu be secured
        :return:
        """
        return value if type(value) is str else ''.join(value)

    @staticmethod
    def plain_mode(value):
        """
        Method will ensure that value will be valid text in latex
        no equation
        :type value: str or list
        :param value: value tu be secured
        :return:
        """
        value = value if type(value) is str else ''.join(value)
        value = value \
            .replace('_', '\\_') \
            .replace("'", "'") \
            .replace('->', '$\\rightarrow$') \
            .replace('<-', '$\\leftarrow$') \
            .replace('\n\n', '\n')
        return str(value)

    @staticmethod
    def name_mode(value):
        """
        Method will ensure that value will be valid name in latex
        This method will replace all characters except a-z A-Z 0-9 and -
        :type value: str or list
        :param value: value tu be secured
        :return:
        """
        value = value if type(value) is str else ''.join(value)
        return re.sub('[^a-zA-Z0-9-]', '', value)