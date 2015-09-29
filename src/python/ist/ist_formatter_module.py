# encoding: utf-8
# author:   Jan Hybs

#
import json, datetime

from ist.formatters.json2html import HTMLFormatter
from ist.formatters.json2latex import LatexFormatter
from ist.nodes import TypedList
from ist.utils.htmltree import htmltree


class ProfilerJSONDecoder(json.JSONDecoder):
    def decode(self, json_string):
        default_obj = super(ProfilerJSONDecoder, self).decode(json_string)
        lst = TypedList()
        lst.parse(default_obj)
        return lst


class ISTFormatter(object):
    """
    Class for formatting json to other formats
    """

    @staticmethod
    def json2latex(input_file='examples/example.json', output_file='../../docs/input_reference_red.tex'):
        """
        Method converts given json file to latex output file
        :param input_file:
        :param output_file:
        :return:
        """
        with open(input_file, 'r') as fp:
            json_object = json.load(fp, encoding="utf-8", cls=ProfilerJSONDecoder)

        latex_result = ''.join(LatexFormatter.format(json_object))
        with open(output_file, 'w') as fp:
            fp.write(latex_result)


    @staticmethod
    def json2html(input_file, output_file, focus_element_id='root', skip_block_creation=[]):
        """
        Method converts given input file to single html output file
        :param input_file:  input json file
        :param output_file: output html file
        :param focus_element_id: id of element which will be visible, default root
        :param skip_block_creation: list of items which won't be created:
         [title, button-control, left-list, ist, right-list]
        :return:
        """
        with open(input_file, 'r') as fp:
            json_object = json.load(fp, encoding="utf-8", cls=ProfilerJSONDecoder)

        html_content = HTMLFormatter.format(json_object)
        html_nav_abc = HTMLFormatter.abc_navigation_bar(json_object)
        html_nav_tree = HTMLFormatter.tree_navigation_bar(json_object)

        # show specified element by given id
        for child in html_content.current():
            if child.attrib['id'].lower() == focus_element_id.lower():
                child.attrib['class'] = child.attrib['class'].replace('hidden', '')
                print "Element {:s} is displayed".format(child.attrib['id'])
                break

        max_cols = 12
        if 'right-list' not in skip_block_creation or 'left-list' not in skip_block_creation:
            max_cols -= 3

        html_body = htmltree('body')
        generated = 'Generated {:s}'.format(datetime.datetime.today().strftime('%d-%m-%Y %X'))
        html_body.span(generated, id='info-generated')
        with html_body.open('div', '', cls='jumbotron', id='top'):
            with html_body.open('div', '', cls='container'):

                if 'title' not in skip_block_creation:
                    with html_body.open('h1', 'Flow123d '):
                        html_body.tag('small', 'input reference')

                if 'search' not in skip_block_creation:
                    with html_body.open('div', cls='form-group has-default has-feedback search-wrapper'):
                        html_body.tag('span', cls='glyphicon glyphicon-search form-control-feedback')
                        html_body.tag('input', type='text', cls='form-control', id='search',
                                      placeholder='type to search')

                with html_body.open('div', cls='row'):
                    if 'right-list' not in skip_block_creation or 'left-list' not in skip_block_creation:
                        with html_body.open('div', cls="col-md-3 tree-list"):
                            with html_body.open('ul', cls="nav nav-tabs", role="tablist"):
                                classes = ['', 'active']

                                if 'left-list' not in skip_block_creation:
                                    with html_body.open('li', cls=classes.pop(), role="presentation"):
                                        html_body.tag('a', 'Tree', attrib={
                                            'role': 'rab',
                                            'aria-controls': 'tree-view',
                                            'data-toggle': 'tab'
                                        })

                                if 'right-list' not in skip_block_creation:
                                    with html_body.open('li', role="presentation", cls=classes.pop()):
                                        html_body.tag('a', 'Abcd', attrib={
                                            'role': 'rab',
                                            'aria-controls': 'abc-view',
                                            'data-toggle': 'tab'
                                        }, )

                            with html_body.open('div', cls='tab-content'):
                                classes = ['tab-pane', 'tab-pane active']
                                if 'left-list' not in skip_block_creation:
                                    with html_body.open('div', role='tabpanel', cls=classes.pop(), id='tree-view'):
                                        html_body.add(html_nav_tree.current())
                                if 'right-list' not in skip_block_creation:
                                    with html_body.open('div', role='tabpanel', cls=classes.pop(), id='abc-view'):
                                        html_body.add(html_nav_abc.current())

                    if 'ist' not in skip_block_creation:
                        with html_body.open('div', cls='col-md-{:d} input-reference'.format(max_cols)):

                            if 'button-control' not in skip_block_creation:
                                with html_body.open('div', id='button-control', cls='row'):
                                    with html_body.open('div', cls='col-md-12'):
                                        with html_body.open('div', id='btn-filter-one-wrapper'):
                                            html_body.tag('input', '', attrib={
                                                'type': 'checkbox',
                                                'class': 'btn btn-default',
                                                'id': 'btn-filter-one',
                                                'data-toggle': 'toggle',
                                                'data-on': 'Single-item',
                                                'data-off': 'Multi-item',
                                                'checked': 'checked'
                                            })
                                        with html_body.open('div', cls='btn-group filter-btns'):
                                            btn_cls = dict()

                                            btn_cls['data-type'] = 'record'
                                            btn_cls['class'] = 'btn btn-warning btn-filter'
                                            html_body.tag('a', 'Records', btn_cls)

                                            btn_cls['data-type'] = 'abstract-record'
                                            btn_cls['class'] = 'btn btn-success btn-filter'
                                            html_body.tag('a', 'Abstract records', btn_cls)

                                            btn_cls['data-type'] = 'selection'
                                            btn_cls['class'] = 'btn btn-info btn-filter'
                                            html_body.tag('a', 'Selections', btn_cls)

                            with html_body.open('div', cls='row'):
                                with html_body.open('a', id='top-link-block', title='Scroll to top',
                                                    href='#input-reference', cls='well well-sm'):
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

        with open(output_file, 'w') as fp:
            fp.write(r'<!DOCTYPE html>')
            fp.write(html.dump())


if __name__ == '__main__':
    formatter = ISTFormatter()
    formatter.json2latex(
        input_file='/home/jan-hybs/Dokumenty/Flow123d-python-utils/src/ist/examples/example.json',
        output_file='/home/jan-hybs/Dokumenty/Smartgit-flow/flow123d/doc/reference_manual/input_reference.tex'
    )

    formatter.json2html(
        input_file='/home/jan-hybs/Dokumenty/Flow123d-python-utils/src/ist/examples/example.json',
        output_file='/home/jan-hybs/Dokumenty/Smartgit-flow/flow123d/doc/reference_manual/input_reference.html',
        focus_element_id='root',
        # skip_block_creation=['title', 'left-list', 'right-list', 'search', 'button-control']
    )
