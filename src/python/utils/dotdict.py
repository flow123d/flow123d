# encoding: utf-8
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

            # name = name.replace('_', '-')
            # if name in self:
            #     return self[name]
            # raise
