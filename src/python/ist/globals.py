# encoding: utf-8
# author:   Jan Hybs
from ist.utils.utils import TypedList


class Globals(object):
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
    def search_in_element(element, value):
        from ist.nodes import Record

        if type(element) is Record:
            return getattr(element, value, None)

        if type(element) is TypedList:
            for item in element:
                try:
                    if item.get('key', 'name').lower() == value.lower():
                        return item
                except:
                    pass
            return None

        # last resort
        return getattr(element, value, None)


    @staticmethod
    def get_url_by_name(label, type=''):
        from ist.utils.htmltree import htmltree
        '''
        constructs and returns tuple (name, type, link) from given name and label
        :param label: name#field where field is optional
        :param type:
        :return:
        '''
        parts = label.split("#")
        name = parts[0]
        fields = parts[1:]

        for (id, item) in Globals.items.iteritems():
            try:
                if item.get('type_name', 'name', 'id').lower() == name.lower():

                    # only link to item?
                    if not fields:
                        return item.get_name(), item.get_type(), '#' + htmltree.chain_values(item.get_name())

                    # link to item's field
                    if fields:
                        curr = item
                        for field in fields:
                            find = Globals.search_in_element(curr, field)
                            if not find:
                                print 'cannot find {} on {}'.format(field, curr)
                                return None
                            # next level
                            curr = find

                        return curr.get_name(), curr.get_type(), '#' + htmltree.chain_values(curr.get_name(), item.get_name())
            except:
                # no such attribute
                pass

        return (None, None, None)