# encoding: utf-8
# author:   Jan Hybs

#
from __future__ import absolute_import

import json, datetime
import re
from ist.formatters.json2html import HTMLFormatter
from ist.formatters.json2latex import LatexFormatter
from ist.globals import Globals
from ist.utils.htmltree import htmltree

from utils.logger import Logger


class ProfilerJSONDecoder(json.JSONDecoder):
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

        latex_result = ''.join(LatexFormatter.format(items))
        with open(output_file, 'w') as fp:
            fp.write(latex_result)


    @staticmethod
    def json2html(items, output_file, focus_element_id='root', skip_block_creation=[],
                  multi_mode=[False, False, True], info=None):
        """
        Method converts given input file to single html output file
        :param items:  list of parsed IST items
        :param output_file: output html file
        :param focus_element_id: id of element which will be visible, default root
        :param skip_block_creation: list of items which won't be created:
         [title, button-control, tree-list, ist, abcd-list]
        :return:
        """
        multi_mode.reverse()
        html_content = HTMLFormatter.format(items)
        html_nav_abc = HTMLFormatter.abc_navigation_bar(items)
        html_nav_tree = HTMLFormatter.tree_navigation_bar(items)
        skip_block_creation.append('tree-list')

        # show specified element by given id
        for child in html_content.current():
            if child.attrib['id'].lower() == focus_element_id.lower():
                child.attrib['class'] = child.attrib['class'].replace('hidden', '')
                break

        max_cols = 12
        if 'abcd-list' not in skip_block_creation or 'tree-list' not in skip_block_creation:
            max_cols -= 3

        html_body = htmltree('body')
        generated = 'Generated {:s}'.format(datetime.datetime.today().strftime('%d-%m-%Y %X'))
        html_body.span(generated, id='info-generated')
        with html_body.open('div', '', cls='jumbotron', id='top'):
            with html_body.open('div', '', cls='container'):

                if 'title' not in skip_block_creation:
                    with html_body.open('h1', 'Flow123d '):
                        version = 'input reference' if not info else info['version']
                        if version.startswith('0.0.'):
                            version = version.replace('0.0.', '')
                        version = re.sub(r'_+', ' ', version)
                        html_body.tag('small', version)

                if 'search' not in skip_block_creation:
                    with html_body.open('div', cls='form-group has-default has-feedback search-wrapper'):
                        html_body.tag('span', cls='glyphicon glyphicon-search form-control-feedback')
                        html_body.tag('input', type='text', cls='form-control', id='search',
                                      placeholder='type to search')

                with html_body.open('div', cls='row'):
                    if 'abcd-list' not in skip_block_creation or 'tree-list' not in skip_block_creation:
                        with html_body.open('div', cls="col-md-3 tree-list"):
                            with html_body.open('ul', cls="nav nav-tabs", role="tablist"):
                                classes = ['', 'active']

                                if 'tree-list' not in skip_block_creation:
                                    with html_body.open('li', cls=classes.pop(), role="presentation"):
                                        html_body.tag('a', 'Tree', attrib={
                                            'role': 'rab',
                                            'aria-controls': 'tree-view',
                                            'data-toggle': 'tab'
                                        })

                                if 'abcd-list' not in skip_block_creation:
                                    with html_body.open('li', role="presentation", cls=classes.pop()):
                                        html_body.tag('a', 'Abc', attrib={
                                            'role': 'rab',
                                            'aria-controls': 'abc-view',
                                            'data-toggle': 'tab'
                                        }, )

                            with html_body.open('div', cls='tab-content'):
                                classes = ['tab-pane', 'tab-pane active']
                                if 'tree-list' not in skip_block_creation:
                                    with html_body.open('div', role='tabpanel', cls=classes.pop(), id='tree-view'):
                                        html_body.add(html_nav_tree.current())
                                if 'abcd-list' not in skip_block_creation:
                                    with html_body.open('div', role='tabpanel', cls=classes.pop(), id='abc-view'):
                                        html_body.add(html_nav_abc.current())

                    if 'ist' not in skip_block_creation:
                        with html_body.open('div', cls='col-md-{:d} input-reference'.format(max_cols)):

                            if 'button-control' not in skip_block_creation:
                                with html_body.open('div', id='button-control', cls='row'):
                                    with html_body.open('div', cls='col-md-12'):
                                        with html_body.open('div', id='btn-filter-one-wrapper'):
                                            single_filter = html_body.tag('input', '', attrib={
                                                'type': 'checkbox',
                                                'class': 'btn btn-default',
                                                'id': 'btn-filter-one',
                                                'data-toggle': 'toggle',
                                                'data-on': 'Single-item',
                                                'data-off': 'Multi-item'
                                            })
                                            if not max(multi_mode):
                                                single_filter.attrib['checked'] = 'checked'

                                        with html_body.open('div', cls='btn-group filter-btns'):
                                            btn_cls = dict()

                                            btn_cls['data-type'] = 'record'
                                            btn_cls['class'] = 'btn btn-warning btn-filter'
                                            btn_cls['class'] += ' active' if multi_mode.pop() else ''
                                            html_body.tag('a', 'Records', btn_cls.copy())

                                            btn_cls['data-type'] = 'abstract-record'
                                            btn_cls['class'] = 'btn btn-success btn-filter'
                                            btn_cls['class'] += ' active' if multi_mode.pop() else ''
                                            html_body.tag('a', 'Abstract records', btn_cls.copy())

                                            btn_cls['data-type'] = 'selection'
                                            btn_cls['class'] = 'btn btn-info btn-filter'
                                            btn_cls['class'] += ' active' if multi_mode.pop() else ''
                                            html_body.tag('a', 'Selections', btn_cls.copy())

                            with html_body.open('div', cls='row'):
                                with html_body.open('a', id='top-link-block', title='Scroll to top',
                                                    href='#top', cls='well well-sm'):
                                    html_body.span(cls='glyphicon glyphicon-menu-up')
                                html_body.add(html_content.current())

        html_head = htmltree('head')

        html_head.tag('title', 'Flow123d input reference')
        html_head.style('css/main.css')
        html_head.style('css/bootstrap.min.css')
        html_head.style('css/bootstrap-toggle.min.css')
        html_head.style('css/katex.min.css')

        html_body.script('js/jquery-2.1.3.min.js')
        html_body.script('js/bootstrap.min.js')
        html_body.script('js/bootstrap-toggle.min.js')
        html_body.script('js/katex.min.js')
        html_body.script('js/main.js')

        html = htmltree('html')
        html.add(html_head.current())
        html.add(html_body.current())

        # try:
        #     Logger.instance().info('Trying module BeautifulSoup')
        #     from BeautifulSoup import BeautifulSoup
        #
        #     Logger.instance().info('Converting to html string')
        #     html_string = html.dump()
        #     Logger.instance().info('Parsing html string')
        #     soup = BeautifulSoup(html_string)
        #     Logger.instance().info('Prettifying html string')
        #     html_pretty = soup.prettify()
        #
        #     Logger.instance().info('Writing to file')
        #     with open(output_file, 'w') as fp:
        #         fp.write(r'<!DOCTYPE html>')
        #         fp.write(html_pretty)
        #     Logger.instance().info('File created')
        # except ImportError as e:
        #     Logger.instance().warning('Import error', exc_info=e)
        Logger.instance().info('Using module ElementTree')
        import xml.etree.ElementTree as ET

        Logger.instance().info('Converting to html string')
        html_string = html.dump()
        Logger.instance().info('Writing to file')
        with open(output_file, 'w') as fp:
            fp.write(r'<!DOCTYPE html>')
            fp.write(ET.tostring(html.root, method='html'))
        Logger.instance().info('File created')
