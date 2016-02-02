# encoding: utf-8
# author:   Jan Hybs


class Field(object):
    """
    Helper class for Single Field in IST node
    """
    def __init__(self, name, value=None, type=str, save_as=None):
        self.name = name
        self.value = value
        self.save_as = save_as or self.name
        self.type = type if value is None else value.__class__

    def __repr__(self):
        if issubclass(self.type, TypedList):
            return "<'{self.name}[{self.type}]'>".format(self=self)
        return "<'{self.name}'>".format(self=self)


class TypedList(list):
    def __init__(self, cls=None, *args, **kwargs):
        self.cls = cls
        list.__init__(self, *args, **kwargs)

    def parse(self, lst):
        """
        method parses given object
        :param lst: object to parse
        :return: self
        """
        from ist.nodes import registered_nodes

        if self.cls is None:
            for o in lst:
                input_type = o.get('input_type', '')
                cls = registered_nodes.get(input_type, None)
                if cls:
                    self.append(cls().parse(o))
                else:
                    self.append(o)
        else:
            for o in lst:
                self.append(self.cls().parse(o))
        return self

    def copy(self):
        """
        Return copy of this instance
        :return:
        """
        return TypedList(self.cls)


class AttributeDict(dict):
    """
    Helper dict attribute class
    """
    def __init__(self, cls=None, *args, **kwargs):
        self.cls = cls
        dict.__init__(self, *args, **kwargs)

    def parse(self, lst):
        """
        method parses given object
        :param lst: object to parse
        :return: self
        """
        for attribute in lst:
            self[attribute] = lst[attribute]

        return self

    def copy(self):
        """
        Return copy of this instance
        :return:
        """
        return AttributeDict(self.cls)