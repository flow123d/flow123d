#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs


class DotDict(dict):
    """ Helper class which enabled field access via dot
    General dict class doesn't allow accessing its keys any
    other way than dict['key_name']

    By overloading __getattribute__ it is possible to use
    dict.key_name

    More, if fields are in somewhat readable format (such as key-name)
    by using dict.key_name will still get desired field
    """

    def __init__(self, *a, **kw):
        dict.__init__(self, *a, **kw)

    def __getattribute__(self, name):
        try:
            return dict.__getattribute__(self, name)
        except AttributeError:
            if name in self:
                return self[name]


class Map(dict):
    """
    Example:
    m = Map({'first_name': 'Eduardo'}, last_name='Pool', age=24, sports=['Soccer'])
    """

    def __init__(self, *args, **kwargs):
        super(Map, self).__init__(*args, **kwargs)
        for arg in args:
            if isinstance(arg, dict):
                for k, v in list(arg.items()):
                    self[k] = v

        if kwargs:
            for k, v in list(kwargs.items()):
                self[k] = v

    def __getattr__(self, attr):
        return self.get(attr)

    def __setattr__(self, key, value):
        self.__setitem__(key, value)

    def __setitem__(self, key, value):
        super(Map, self).__setitem__(key, value)
        self.__dict__.update({key: value})

    def __delattr__(self, item):
        self.__delitem__(item)

    def __delitem__(self, key):
        super(Map, self).__delitem__(key)
        del self.__dict__[key]
