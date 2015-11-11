# encoding: utf-8
# author:   Jan Hybs


'''
MdLatex Extension for Python-Markdown
======================================

Converts [[type_value]] to relative links.
'''



import cgi
from markdown.extensions import Extension
from markdown.inlinepatterns import Pattern, SimpleTagPattern, BacktickPattern
from markdown.util import etree
import re


class MdLatexExtension (Extension):
    def __init__ (self, *args, **kwargs):
        self.config = {
        }

        super (MdLatexExtension, self).__init__ (*args, **kwargs)

    def extendMarkdown (self, md, md_globals):
        self.md = md

        # append to end of inline patterns
        MD_LATEX_RE = r'\(\((.*?)\)\)'  # match anything in (( ))
        md_latex_pattern = MdLatex (MD_LATEX_RE, self.getConfigs ())
        md_latex_pattern.md = md
        md.inlinePatterns.add ('mdlatex', md_latex_pattern, "_begin")

        # md.treeprocessors.add(
        # "footnote", FootnoteTreeprocessor(), "_begin"
        # )


class MdLatex (Pattern):
    def __init__ (self, pattern, config):
        super (MdLatex, self).__init__ (pattern)
        self.config = config

    def handleMatch (self, m):
        if m.group (2).strip ():
            latex_expression = m.group (2).strip ()
            span = etree.Element ('span')
            span.text = latex_expression
            span.set ('class', 'md-expression')

            return span
        else:
            return ''


def makeExtension (*args, **kwargs):
    return MdLatexExtension (*args, **kwargs)


# from markdown.treeprocessors import Treeprocessor
#
# class FootnoteTreeprocessor(Treeprocessor):
# def run(self, root):
# root.text = '!!!!!!!!!!!!!!!'


class MdLatexSupport (object):
    def __init__ (self):
        self.latex = []

    def _match_prepare (self, m):
        self.latex.append (m.group (1))
        return '(({}))'.format (len (self.latex) - 1)

    def _match_finish (self, m):
        index = int (m.group (1))
        latex = cgi.escape(self.latex[index])
        return '<span class="md-expression">{{{}}}</span>'.format (latex)

    def prepare (self, html):
        regex = r'\(\(([\w\d\s{}\[\]\\\?;$^?_/+!&=*<>~ ()-]*?)\)\)'

        secured_html = re.sub (regex, self._match_prepare, html, re.S | re.M | re.UNICODE)
        return secured_html

    def finish (self, html):
        regex = r'\(\(([0-9]+)\)\)'

        secured_html = re.sub (regex, self._match_finish, html, re.S | re.M | re.UNICODE)
        return secured_html