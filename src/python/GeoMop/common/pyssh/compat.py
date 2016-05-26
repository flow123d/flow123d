# -*- coding: utf-8 -*-

import six

text_type = six.text_type
binary_type = six.binary_type

if six.PY2:
    def to_bytes(data, encoding="utf-8"):
        if isinstance(data, text_type):
            return data.encode(encoding)
        return data

else:
    def to_bytes(data, encoding="utf-8"):
        if isinstance(data, text_type):
            return bytes(data, encoding)
        return data

def to_text(data, encoding="utf-8"):
    if isinstance(data, binary_type):
        return data.decode(encoding)
    return data
