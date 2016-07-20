#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs

from __future__ import absolute_import
from utils.logger import Logger


class Globals(object):
    """
    Global class object which stores references and all objects on memory for later use
    :type items             : dict[ist.base.Parsable]
    """
    items = { }
    names = {
        'record': 'type_name',
        'r': 'type_name',
        'abstractrecord': 'name',
        'a': 'name',
        'ar': 'name',
        'selection': 'name',
        's': 'name',
        '': 'name'
    }

    @staticmethod
    def iterate(type=''):
        from ist.base import InputType
        """
        :rtype : list[ist.base.Parsable]
        """
        if not type:
            return Globals.items.values()
        elif type.lower() in ('r', 'record'):
            return [x for x in Globals.items.itervalues() if getattr(x, 'input_type', InputType.UNKNOWN) == InputType.RECORD]
        elif type.lower() in ('s', 'selection'):
            return [x for x in Globals.items.itervalues() if getattr(x, 'input_type', InputType.UNKNOWN) == InputType.SELECTION]
        elif type.lower() in ('a', 'ar', 'abstract'):
            return [x for x in Globals.items.itervalues() if getattr(x, 'input_type', InputType.UNKNOWN) == InputType.ABSTRACT_RECORD]

    @staticmethod
    def get_url_by_name(label, item_type=''):
        """
        constructs and returns tuple (name, type, link) from given name and label
        :param label: name#field where field is optional
        :param item_type:
        :rtype : list[ist.base.Parsable]
        """
        parts = label.split("#")
        name = parts[0]
        field = parts[1] if len(parts) > 1 else None

        for item in Globals.iterate(item_type):
            possibilities = item.gets('name', 'id')
            if getattr(item, 'attributes', None):
                possibilities.extend(item.attributes.gets('link_name'))
            possibilities = [str(p).lower() for p in possibilities if p]

            if name.lower() in possibilities:
                # only link to item?
                if not field:
                    return item, None
                else:
                    return item, item.get_fields(field)

        return None, None, None

    @staticmethod
    def save(key, item):
        if key in Globals.items:
            Logger.instance().info('duplicate key %s' % key)
            for i in range(2, 100):
                new_key = '{}-{}'.format(key, i)
                if new_key not in Globals.items:
                    Logger.instance().info('For key %s assigned %s' % (key, new_key))
                    Globals.items[new_key] = item
                    return new_key
        Globals.items[key] = item
        return key


class FormatMode(object):
    LATEX_MODE = 1
    HTML_MODE = 2
    format_mode = 0
