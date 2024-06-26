#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

"""
MDStrike Extension for Python-Markdown
======================================

Converts ~~expression~~ to line-through text <del /> tag is used
"""

from markdown.extensions import Extension
from markdown.inlinepatterns import Pattern
import xml.etree.ElementTree as etree


class StrikeThroughExtension(Extension):
    """
    Class StrikeThroughExtension is module for md allowing to strike through the text
    """

    def __init__(self, *args, **kwargs):
        self.md = None
        self.config = {
        }

        super(StrikeThroughExtension, self).__init__(*args, **kwargs)

    def extendMarkdown(self, md):
        self.md = md

        # append to end of inline patterns
        # WIKILINK_RE = r'\[\[([\w0-9_ -]+)\]\]'
        WIKILINK_RE = r'~~([\w0-9_#-]+)~~'
        wikilinkPattern = StrikeThroughPattern(WIKILINK_RE, {})
        wikilinkPattern.md = md
        md.inlinePatterns.register(wikilinkPattern, 'mdstrikethrough', 175)


class StrikeThroughPattern(Pattern):
    """
    Class StrikeThroughPattern handles single match in md
    """

    def __init__(self, pattern, config):
        super(StrikeThroughPattern, self).__init__(pattern)
        self.config = config

    def handleMatch(self, m):
        label = m.group(2).strip()
        if label:
            span = etree.Element('del')
            span.text = label
            return span
        return ''


def makeExtension(*args, **kwargs):
    """
    register extension
    :param args:
    :param kwargs:
    :return:
    """
    return StrikeThroughExtension(*args, **kwargs)
