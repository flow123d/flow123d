#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

from __future__ import absolute_import
import os
import json
import time
from ist.formatters.extensions.json2latex2 import LatexFormatter
from ist.formatters.json2html import HTMLFormatter
from ist.globals import FormatMode
from ist.utils.htmltree import htmltree
from utils.logger import Logger


class ProfilerJSONDecoder(json.JSONDecoder):
    """
    Class ProfilerJSONDecoder is custom json decoder
    """

    def decode(self, json_string):
        default_obj = super(ProfilerJSONDecoder, self).decode(json_string)
        return default_obj


class ISTFormatter(object):
    """
    Class for formatting json to other formats
    """

    @staticmethod
    def json2latex(items, output_file='../../docs/input_reference_red.tex', info=None):
        """
        Method converts given json file to latex output file
        :param items:  list of parsed IST items
        :param output_file:
        :return:
        """
        FormatMode.format_mode = FormatMode.LATEX_MODE
        latex_result = LatexFormatter.format(items).to_string()
        with open(output_file, 'w') as fp:
            fp.write(latex_result)

    @staticmethod
    def tree2html(element, try_pretty=0):
        import xml.etree.ElementTree as ET
        html_string = ET.tostring(element.root, method='html')
        if not try_pretty:
            return html_string

        try:
            from BeautifulSoup import BeautifulSoup
            soup = BeautifulSoup(html_string)
            Logger.instance().info('Prettifying html string')
            html_pretty = soup.prettify()
            return html_pretty
        except:
            pass

        return html_string

    @staticmethod
    def json2html(items, output_file, focus_element_id='root', skip_block_creation=[],
                  multi_mode=[1, 0, 0, 0], info=None):
        """
        Method converts given input file to single html output file
        :param items:  list of parsed IST items
        :param output_file: output html file
        :param focus_element_id: id of element which will be visible, default root
        :param skip_block_creation: list of items which won't be created:
         [title, button-control, tree-list, ist, abcd-list]
        :return:
        """
        FormatMode.format_mode = FormatMode.HTML_MODE
        template_file = os.path.join(os.path.dirname(__file__), 'formatters/templates/template.html')
        template_string = open(template_file, 'r').read()

        html_items = HTMLFormatter.format(items)
        html_items_string = ISTFormatter.tree2html(html_items)
        # html_items_string = "foo"

        html_navigation = HTMLFormatter.abc_navigation_bar(items)
        html_navigation_string = ISTFormatter.tree2html(html_navigation)

        template_string = template_string.replace("@GENERATED@", time.strftime("%d-%m-%Y %H:%M:%S"))
        template_string = template_string.replace("@IST@", html_items_string)
        template_string = template_string.replace("@NAVIGATION@", html_navigation_string)

        with open(output_file, 'w') as fp:
            fp.write(template_string)
