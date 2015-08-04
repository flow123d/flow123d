# encoding: utf-8
# author:   Jan Hybs

#
import json
from ist.formatters.extensions.md_latex import MdLatexSupport

from ist.formatters.json2html import HTMLFormatter
from ist.formatters.json2latex import LatexFormatter
from ist.formatters.markdown2html import markdown2html
from ist.nodes import TypedList
from ist.utils.htmltree import htmltree


class ProfilerJSONDecoder(json.JSONDecoder):
    def decode(self, json_string):
        default_obj = super(ProfilerJSONDecoder, self).decode(json_string)
        lst = TypedList()
        lst.parse(default_obj)
        return lst

class ISTFormatter (object):
    """
    Class for formatting json to other formats
    """
    def json2latex(self, input_file='examples/example.json', output_file='../../docs/input_reference_red.tex'):
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


    def json2html(self, input_file='examples/example.json', output_file='../../docs/index.html'):
        """
        Method converts given input file to single html output file
        :param input_file:
        :param output_file:
        :return:
        """
        with open(input_file, 'r') as fp:
            json_object = json.load(fp, encoding="utf-8", cls=ProfilerJSONDecoder)

        html_content = HTMLFormatter.format(json_object)
        html_nav_abc = HTMLFormatter.abc_navigation_bar(json_object)
        html_nav_tree = HTMLFormatter.tree_navigation_bar(json_object)

        html_body = htmltree('body')
        with html_body.open('div', '', { 'class': 'jumbotron' }):
            with html_body.open('div', '', { 'class': 'container' }):
                with html_body.open('h1', 'Flow123d '):
                    html_body.tag('small', 'input reference')

                with html_body.open('div', attrib={ 'class': 'btn-group', 'data-toggle': 'buttons' }):
                    btn_cls = dict({ 'class': 'btn btn-default btn-filter' })
                    btn_cls['data-type'] = 'record'
                    html_body.tag('a', 'Records', btn_cls)
                    btn_cls['data-type'] = 'abstract-record'
                    html_body.tag('a', 'Abstract records', btn_cls)
                    btn_cls['data-type'] = 'selection'
                    html_body.tag('a', 'Selections', btn_cls)

                html_body.tag('a', 'Filter one', attrib={ 'class': 'btn btn-default btn-filter-one' })

                with html_body.open('div', attrib={ 'class': 'col-md-2 tree-list' }):
                    html_body.add(html_nav_tree.current())
                with html_body.open('div', attrib={ 'class': 'col-md-8 input-reference' }):
                    html_body.add(html_content.current())
                    first_record = html_body.current()._children[0]._children[0]
                    first_record.attrib['class'] = first_record.attrib['class'].replace('hidden', '')
                with html_body.open('div', attrib={ 'class': 'col-md-2 abc-list' }):
                    html_body.add(html_nav_abc.current())

        html_head = htmltree('head')

        html_head.tag('title', 'Flow123d input reference')
        html_head.style('css/main.css')
        html_head.style('css/bootstrap.min.css')
        html_head.style('css/katex.min.css')

        html_body.script('js/jquery-2.1.3.min.js')
        html_body.script('js/bootstrap.min.js')
        html_body.script('js/katex.min.js')
        html_body.script('js/main.js')

        html = htmltree('html')
        html.add(html_head.current())
        html.add(html_body.current())

        with open(output_file, 'w') as fp:
            fp.write(r'<!DOCTYPE html>')
            fp.write(html.dump())


if __name__ == '__main__':
    print 'converting'
    formatter = ISTFormatter ()
    formatter.json2latex(
        input_file='/home/jan-hybs/Dokumenty/Smartgit-flow/flow123d/doc/reference_manual/input_reference.json',
        output_file='/home/jan-hybs/Dokumenty/Smartgit-flow/flow123d/doc/reference_manual/input_reference.tex'
    )
