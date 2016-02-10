# encoding: utf-8
# author:   Jan Hybs


'''
MdLatex Extension for Python-Markdown
======================================

Converts [[type_value]] to relative links.
'''

from __future__ import absolute_import
from __future__ import unicode_literals
from markdown.extensions import Extension
from markdown.inlinepatterns import Pattern
from markdown.util import etree
import re
from ist.globals import Globals
from utils.logger import Logger


def build_url(label, base, end):
    """ Build a url from the label, a base, and an end. """
    clean_label = re.sub(r'([ ]+_)|(_[ ]+)|([ ]+)', '_', label)
    return '%s%s%s' % (base, clean_label, end)


class MdLinkExtension(Extension):
    def __init__(self, *args, **kwargs):
        self.md = None
        self.config = {
        }

        super(MdLinkExtension, self).__init__(*args, **kwargs)

    def extendMarkdown(self, md, md_globals):
        self.md = md

        # append to end of inline patterns
        # WIKILINK_RE = r'\[\[([\w0-9_ -]+)\]\]'
        WIKILINK_RE = r'\[\[([\w0-9_#:-]+)\]\]'
        wikilinkPattern = MdLinks(WIKILINK_RE, { })
        wikilinkPattern.md = md
        md.inlinePatterns.add('mdlinks', wikilinkPattern, "<not_strong")


class MdLinks(Pattern):
    def __init__(self, pattern, config):
        super(MdLinks, self).__init__(pattern)
        self.config = config

    def handleMatch(self, match):
        if match.group(2).strip():
            label = match.group(2).strip()
            link_type = ''
            pos = label.find('#')
            if pos != -1:
                opts = ('r', 'record', 's', 'selection', 'a', 'abstract', 'ar')
                if label[:pos] in opts:
                    link_type = label[:pos]
                    label = label[pos+1:]

            element = self.build_element(link_type, label)
            return element
        else:
            return ''

    @staticmethod
    def build_element(link_type, label):

        # find item which is desired
        link_text = None
        if label.find(':') != -1:
            link_text = label[label.find(':') + 1:]
            label = label[:label.find(':')]

        result = Globals.get_url_by_name(label, link_type)
        item = result[0]
        item_field = result[1]
        a = etree.Element('a')

        if item_field:
            a.text = link_text or item_field.href_name
            a.set('href', '#{item_field.href_id}'.format(item_field=item_field))

        elif item:
            a.text = link_text or item.href_name
            a.set('href', '#{item.href_id}'.format(item=item))
        else:
            Logger.instance().warning('Link not found %s %s' % (link_type, label) )
            return ''
        return a


def makeExtension(*args, **kwargs):
    return MdLinkExtension(*args, **kwargs)
