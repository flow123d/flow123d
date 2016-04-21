#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

from __future__ import absolute_import
import xml.etree.ElementTree as ET


class Html2Latex(object):
    """
    Class which based on given element (html/string) can produce latex format
    """
    list_types = {
        'ul': 'itemize',
        'ol': 'enumerate'
    }

    def __init__(self, element):
        if type(element) in (str, unicode):
            tree = ET.fromstring('<html_example>' + element + "</html_example>")
            self.el = tree
        else:
            self.el = element

        from ist.utils.texlist2 import TexList

        self.tex = TexList()

    def extend_children(self):
        """
        Recursive children call
        """
        for child in self.el:
            self.tex.extend(Html2Latex(child).to_latex())

    def tag_is(self, *tags):
        """
        whether current tag is in given tags
        """
        return self.el.tag in tags

    def text(self):
        """ return current text """
        return self.el.text if self.el.text else ''

    def tail(self):
        """ return current tail """
        return self.el.tail if self.el.tail else ''

    def add_text(self):
        if str(self.text).strip():
            with self.tex:
                self.tex.append(self.text())

    def add_tail(self):
        """
        Adds current tail if exists
        """
        if self.tail():
            if self.tail() == '\n':
                self.tex.append(self.tail())
            else:
                with self.tex:
                    self.tex.append(self.tail())

    def get_list_type(self):
        """ helper method for getting list type"""
        return self.list_types.get(self.el.tag, 'itemize')

    #
    # def tag (self):
    # return self.el.tag

    def to_latex(self):
        """
        Method converts this object to latex with recursive calls
        """

        if self.tag_is('p'):
            with self.tex:
                with self.tex:
                    self.tex.append(self.text())
                self.extend_children()
                self.add_tail()
                # self.tex.newline ()

        elif self.tag_is('br'):
            self.tex.append(self.text())
            self.tex.append(r'\\')
            self.extend_children()
            self.add_tail()

        elif self.tag_is('h1'):
            self.tex.append('\\section')
            with self.tex:
                self.tex.append(self.text())
                self.extend_children()
            self.add_tail()

        elif self.tag_is('a'):
            self.tex.extend(LatexHref(self.el).to_latex())

        elif self.tag_is('em'):
            self.tex.append('\\textit')
            with self.tex:
                self.tex.append(self.text())
                self.extend_children()
            self.add_tail()

        elif self.tag_is('strong'):
            self.tex.append('\\textbf')
            with self.tex:
                self.tex.append(self.text())
                self.extend_children()
            self.add_tail()

        elif self.tag_is('ul', 'ol'):
            self.tex.begin(self.get_list_type())
            self.tex.append(self.text())
            self.extend_children()
            self.tex.end(self.get_list_type())
            self.add_tail()

        elif self.tag_is('li'):
            # with self.tex:
            self.tex.append('\\item ')
            self.add_text()
            self.extend_children()
            self.add_tail()


        # so far, code tag will be monospaced only
        elif self.tag_is('code'):
            # self.tex.open_element ('lstlisting')
            # self.tex.newline ()
            self.tex.slash('ttfamily ')
            self.tex.append(self.text().replace('\$', '\$'))
            # self.extend_children ()
            # self.tex.newline ()
            # self.tex.close_element ('lstlisting')
            # self.add_tail ()

        elif self.tag_is('span'):
            self.tex.append(self.text())
            self.extend_children()
            self.add_tail()

        else:
            with self.tex:
                self.tex.append(self.text())
                self.extend_children()
                self.add_tail()

        return self.tex


class LatexHref(Html2Latex):
    """
    Class for latex hrefs
    """

    def to_latex(self):

        # Alink href?
        if self.el.attrib.get('data-href') == 'Alink':
            url = self.el.attrib.get('href')
            if url.startswith('#'):
                url = url[1:]
            text = self.el.attrib.get('text')
            self.tex.macro_alink_(url=url, text=text)
            self.add_tail()
        else:
            # other href
            self.tex.append('\\href')
            self.tex.add(self.el.attrib.get('href'))
            with self.tex:
                self.tex.add(self.text())
                self.extend_children()
            self.add_tail()
        return self.tex
