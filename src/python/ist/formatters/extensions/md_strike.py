#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

"""
MDStrike Extension for Python-Markdown
======================================

Converts ~~expression~~ to line-through text <del /> tag is used
"""

from __future__ import absolute_import
from __future__ import unicode_literals
from markdown.extensions import Extension
from markdown.inlinepatterns import Pattern
from markdown.util import etree


class StrikeThroughExtension(Extension):
    def __init__(self, *args, **kwargs):
        self.md = None
        self.config = {
        }

        super(StrikeThroughExtension, self).__init__(*args, **kwargs)

    def extendMarkdown(self, md, md_globals):
        self.md = md

        # append to end of inline patterns
        # WIKILINK_RE = r'\[\[([\w0-9_ -]+)\]\]'
        WIKILINK_RE = r'~~([\w0-9_#-]+)~~'
        wikilinkPattern = StrikeThroughPattern(WIKILINK_RE, { })
        wikilinkPattern.md = md
        md.inlinePatterns.add('mdstrikethrough', wikilinkPattern, "<not_strong")


class StrikeThroughPattern(Pattern):
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
    return StrikeThroughExtension(*args, **kwargs)
