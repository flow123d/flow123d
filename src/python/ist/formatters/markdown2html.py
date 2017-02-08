#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

from __future__ import absolute_import

import cgi
import re

import markdown
from xml.etree import ElementTree as ET


class markdown2html(object):
    """
    Class markdown2html is simple helper class for parsing md to html
    """

    def __init__(self):
        self._md_latex = ExpressionPlaceholder()

    def _replace_math(self, md_text):
        """
        Replace math formulas with unique placeholders
        :type md_text: str
        """
        return self._md_latex.prepare(md_text)

    def _replace_placeholders(self, md_text):
        """
        Replace placeholder with original math formulas
        :type md_text: str
        """
        return self._md_latex.finish(md_text)

    def parse2latex(self, md_text):
        from ist.utils.texlist2 import TexList
        reduce_tag = 'div'

        secured_markdown = self._replace_math(md_text)
        latex_secured = TexList.prepare_plain(secured_markdown)

        # apply markdown
        html_secured = markdown.markdown(latex_secured, extensions=[
            'markdown.extensions.sane_lists',
            'markdown.extensions.nl2br',
            'ist.formatters.extensions.md_links',
            'ist.formatters.extensions.md_strike'])
        html_secured = TexList.finish_plain(html_secured)
        html = self._md_latex.finish(html_secured)

        try:
            return ET.fromstring(html)
        except Exception as e:
            return ET.fromstring('<' + reduce_tag + '>' + html + '</' + reduce_tag + '>')

    def parse(self, md_text, reduce_to_tree=False, reduce_tag='div'):
        secured_markdown = self._md_latex.prepare(md_text)
        secured_markdown = cgi.escape(secured_markdown)

        html_secured = markdown.markdown(secured_markdown, extensions=[
            'markdown.extensions.sane_lists',
            'markdown.extensions.nl2br',
            'ist.formatters.extensions.md_links',
            'ist.formatters.extensions.md_strike'])
        html = self._md_latex.finish(html_secured)

        if not reduce_to_tree:
            return html
        try:
            return ET.fromstring(html)
        except Exception as e:
            return ET.fromstring('<' + reduce_tag + '>' + html + '</' + reduce_tag + '>')


class ExpressionPlaceholder(object):
    """
    Class ExpressionPlaceholder is helper class handling replacement in latex
    """

    regex_math_formula = re.compile(r'\(\(([\w\d\s{}\[\]\\\?;$^?_/+!&=*<>~ (),\.-]*?)\)\)')
    regex_placeholder = re.compile(r'\(\(([0-9]+)\)\)')
    result_placeholder = '<span class="md-expression">{{{}}}</span>'

    def __init__(self,
                 match_regex=regex_math_formula,
                 placeholder_regex=regex_placeholder,
                 placeholder_result=result_placeholder):
        self.latex = []
        self.match_regex = match_regex
        self.placeholder_result = placeholder_result
        self.placeholder_regex = placeholder_regex

    def _match_prepare(self, m):
        self.latex.append(m.group(1))
        return '(({}))'.format(len(self.latex) - 1)

    def _match_finish(self, m):
        index = int(m.group(1))
        latex = cgi.escape(self.latex[index])
        return self.placeholder_result.format(latex)

    def prepare(self, html):
        secured_html = self.match_regex.sub(self._match_prepare, html, re.S | re.M | re.UNICODE)
        return secured_html

    def finish(self, html):
        secured_html = self.placeholder_regex.sub(self._match_finish, html, re.S | re.M | re.UNICODE)
        return secured_html
