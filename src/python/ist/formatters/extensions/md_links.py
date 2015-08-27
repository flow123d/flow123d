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


def build_url (label, base, end):
    """ Build a url from the label, a base, and an end. """
    clean_label = re.sub (r'([ ]+_)|(_[ ]+)|([ ]+)', '_', label)
    return '%s%s%s' % (base, clean_label, end)


class MdLinkExtension (Extension):
    def __init__ (self, *args, **kwargs):
        self.config = {
        }

        super (MdLinkExtension, self).__init__ (*args, **kwargs)

    def extendMarkdown (self, md, md_globals):
        self.md = md

        # append to end of inline patterns
        # WIKILINK_RE = r'\[\[([\w0-9_ -]+)\]\]'
        WIKILINK_RE = r'\[\[([\w0-9-]+_[\w0-9_#-]+)\]\]'
        wikilinkPattern = MdLinks (WIKILINK_RE, self.getConfigs ())
        wikilinkPattern.md = md
        md.inlinePatterns.add ('mdlinks', wikilinkPattern, "<not_strong")


class MdLinks (Pattern):
    def __init__ (self, pattern, config):
        super (MdLinks, self).__init__ (pattern)
        self.config = config

    def handleMatch (self, m):
        if m.group (2).strip ():
            label = m.group (2).strip ()
            (type, label) = label.split ("_", 1)
            element = self.build_element (type, label)
            return element
        else:
            return ''

    def build_element (self, type, label):
        if type.lower () == 'attribute':
            p = etree.Element ('p')
            p.text = 'attribute value here'
            return p

        if type.lower () in ('record', 'abstractrecord', 'selection', 'r', 'a', 's', 'ar'):

            # find item which is desired
            (name, type, link) = Globals.get_url_by_name (label)
            if not link:
                return None

            a = etree.Element ('a')
            a.text = '{name} ({type})'.format (name=name, type=type)
            a.set ('data-href', 'Alink')
            a.set ('href', link)
            a.set ('data-ns', 'IT::')
            return a

        print 'unknown type'
        return None


def makeExtension (*args, **kwargs):
    return MdLinkExtension (*args, **kwargs)
