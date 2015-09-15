# encoding: utf-8
# author:   Jan Hybs
from xml.etree import ElementTree as ET
import markdown
from ist.formatters.extensions.md_latex import MdLatexSupport


class markdown2html(object):
    def __init__(self):
        self._md_latex = MdLatexSupport()

    def parse(self, md_text, reduce_to_tree=False, reduce_tag='div'):
        secured_markdown = self._md_latex.prepare(md_text)
        html_secured = markdown.markdown(secured_markdown, extensions=[
            'markdown.extensions.sane_lists',
            'markdown.extensions.nl2br',
            'ist.formatters.extensions.md_links'])
        html = self._md_latex.finish(html_secured)

        if not reduce_to_tree:
            return html
        try:
            return ET.fromstring(html)
        except Exception as e:
            return ET.fromstring('<' + reduce_tag + '>' + html + '</' + reduce_tag + '>')
