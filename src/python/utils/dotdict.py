# encoding: utf-8
# author:   Jan Hybs

class DotDict(dict):
    def __init__(self, *a, **kw):
        dict.__init__(self, *a, **kw)

    def __getattribute__(self, name):
        try:
            return dict.__getattribute__(self, name)
        except AttributeError:
            if name in self:
                return self[name]

            name = name.replace('_', '-')
            if name in self:
                return self[name]
            raise
