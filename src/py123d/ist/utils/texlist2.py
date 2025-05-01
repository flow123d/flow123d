#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
import re
from py123d.ist.formatters.html2latex import Html2Latex
from py123d.ist.formatters.markdown2html import markdown2html


class TexList(list):
    """
    Helper class for creating latex document
    """
    special_chars = [
        ['''\\''', r'{\textbackslash}'],
        [r'{', r'{\{}'],
        [r'}', r'{\}}'],
        [r'%', r'{\%}'],
        [r'$', r'{\$}'],
        [r'#', r'{\#}'],
        [r'&', r'{\&}'],
        [r'_', r'{\_}'],
        [r'^', r'{\^{}}'],
        [r'|', r'{\textbar}'],
        [r'>', r'{\textgreater}'],
        [r'<', r'{\textless}'],
        [r'->', r'{\rightarrow}'],
        [r'<-', r'{\leftarrow}'],
    ]
    m2h = markdown2html()
    _OPEN = '{'
    _CLOSE = '}'
    _SLASH = '\\'
    _IT = 'it'
    _SPACE = ' '
    _BEGIN = 'begin'
    _END = 'end'
    _NEWLINE = '\n'
    _TAB = '\t'
    _TEXTLANGLE = 'textlangle'
    _TEXTRANGLE = 'textrangle'
    """Will create well-formatted latex file but result will be INVALID, for debug purposes only"""
    PRETTY_FORMAT = False

    def __init__(self):
        super(TexList, self).__init__()
        self.levels = {}
        self.l = 0
        self.open_name = ''
        self.append('')


    def add(self, value, t=None):
        if t is None:
            self.append(self._OPEN + str(value) + self._CLOSE)
        else:
            self.append(self._OPEN + str(t(value)) + self._CLOSE)

    def comment(self, value):
        if not self.PRETTY_FORMAT:
            return
        self.append("% " + str(value))

    def to_string(self, fix_newlines=True):
        if not fix_newlines:
            return ''.join(self)

        tmp_list = list()
        nl = False
        for x in self:
            if x == '':
                continue
            if x == self._NEWLINE:
                if nl:
                    tmp_list.append(x)
                nl = False
            else:
                tmp_list.append(x)
                if str(x).strip():
                    nl = True
        return ''.join(tmp_list)

    def _function_call(self, args, mode=None, func=None):
        if func:
            self.slash(func)

        # secure arguments
        args = args if type(args) is list else [args]
        if mode is None:
            mode = [None] * len(args)
        mode = mode if type(mode) is list else [mode]

        for i in range(len(args)):
            self.add(str(args[i]), mode[i])

    def slash(self, value=''):
        """
        Add \value
        """
        if value:
            self.append(self._SLASH + str(value))
        else:
            self.append(self._SLASH)

    def macro_alink(self, item, text=None):
        """
        Adds \Alink{item.href_id}{item.href_name}, will create href
        :type item: ist.nodes.TypeSelection or ist.nodes.TypeRecord or ist.nodes.TypeAbstract or Unicode
        """
        self._function_call(
            func='TypeLink',
            args=[item.href_id, text or item.href_name],
            mode=[self.TYPE_NONE, self.TYPE_PLAIN]
        )

    '''
    def macro_alink_(self, url, text):
        t = TexList()
        t._function_call(
            func='Alink',
            args=[url, text or text],
            mode=[self.TYPE_NAME, self.TYPE_PLAIN]
        )
        self.append(str(t))
    '''

    def macro_hyper_b(self, item, text=None):
        """
        Adds \hyperB{item.href_id}{item.href_name}, will register href
        :type item: Parsable or ist.nodes.TypeSelection or ist.nodes.TypeRecord or ist.nodes.TypeAbstract
        """
        self._function_call(
            func='hyperB',
            args=[item.href_id, text or item.href_name],
            mode=[self.TYPE_NONE, self.TYPE_PLAIN]
        )

    '''
    def macro_add_doc(self, item):
        """
        Adds \hyperB{item.href_id}{item.href_name}, will register href
        :type item: Parsable or ist.nodes.TypeSelection or ist.nodes.TypeRecord or ist.nodes.TypeAbstract
        """
        self._function_call(
            func='AddDoc',
            args=[item.href_name],
            mode=[self.TYPE_PLAIN]
        )
    '''

    def macro_text_lr_angle(self, value, mode=None, italic=True):
        self.slash(self._TEXTLANGLE)
        with self:
            self.append(self._SPACE)
            self._function_call(
                func=self._IT if italic else None,
                args=value,
                mode=mode if mode else self.TYPE_PLAIN
            )
            self.append(self._SPACE)
        self.slash(self._TEXTRANGLE)

    def begin(self, name):
        self._newline()
        self.slash(self._BEGIN)
        self.add(name)
        self._newline()

    def end(self, name):
        self._newline()
        self.slash(self._END)
        self.add(name)
        self._newline()

    def __enter__(self):
        """
        Enter the runtime context related to this object.
        """
        if self.open_name:
            self.levels[self.l] = self.open_name
            self.l += 1
            self.slash(self.open_name)
            self.open_name = ''
        else:
            self.levels[self.l] = False
            self.l += 1
            self.append(self._OPEN)
        return self

    def __exit__(self, exception_type, exception_value, tb):
        """
        Exit the runtime context related to this object.
        """
        if exception_type:
            return False

        self.l -= 1
        if not self.levels[self.l]:
            self.append(self._CLOSE)
        return self

    def _newline(self):
        if not self.PRETTY_FORMAT:
            return
        if self[-1] != self._NEWLINE:
            self.append(self._NEWLINE)

    def _tab(self, n=1):
        if not self.PRETTY_FORMAT:
            return
        self.append(self._TAB * n)

    @classmethod
    def description(cls, value):
        if not value.strip():
            return ''

        html = TexList.m2h.parse2latex(str(value))
        latex = Html2Latex(html)
        result = latex.to_latex()
        str_repr = "".join(result)
        # Fix sentences without space
        str_repr=re.sub(r'([a-z])\.([A-Z])', r'\1. \2', str_repr)
        # split sentences to more lines
        str_repr=re.sub(r'([a-z])\. ([A-Z])', r'\1.\n\2', str_repr)
        return str_repr

    @classmethod
    def equation_mode(cls, value):
        """
        Method will ensure that value will be valid equation in latex
        :type value: str or list
        :param value: value tu be secured
        :return:
        """
        return value if type(value) is str else ''.join(value)

    @classmethod
    def plain_mode(cls, value):
        """
        Method will ensure that value will be valid text in latex
        no equation
        :type value: str or list
        :param value: value tu be secured
        :return:
        """
        value = value if type(value) is str else ''.join(value)

        value = cls.prepare_plain(value)
        return cls.finish_plain(value)

    @classmethod
    def prepare_plain(cls, value):
        # replace all characters with placeholders
        # since some characters contains other non allowed chars
        for i in range(len(cls.special_chars)):
            value = value.replace(cls.special_chars[i][0], '(-(-{}-)-)'.format(i))
        return value

    @classmethod
    def finish_plain(cls, value):
        # replace placeholders with escaped characters
        for i in range(len(cls.special_chars)):
            value = value.replace('(-(-{}-)-)'.format(i), cls.special_chars[i][1])
        return str(value)

    @classmethod
    def name_mode(cls, value):
        """
        Method will ensure that value will be valid name in latex
        This method will replace all characters except a-z A-Z 0-9 and -
        :type value: str or list
        :param value: value tu be secured
        :return:
        """
        value = value if type(value) is str else ''.join(value)
        return re.sub('[^a-zA-Z0-9-]+', '-', value)

    @classmethod
    def auto_mode(cls, values):
        result = list()
        for value in values:
            if value.startswith('{$') and value.endswith('$}'):
                result.append(TexList.TYPE_EQ(value))
            elif value.startswith('\Alink'):
                result.append(value)
            else:
                result.append(TexList.TYPE_PLAIN(value))
        return ''.join(result).rstrip('\\')

    @classmethod
    def none_mode(cls, values):
        if type(values) is str:
            return values
        return ''.join(values)

    def __str__(self):
        return ''.join(self)

    TYPE_NAME = name_mode
    TYPE_EQ = equation_mode
    TYPE_PLAIN = plain_mode
    TYPE_AUTO = auto_mode
    TYPE_NONE = none_mode

    def item_open(self, name):
        self.open_name = name
        return self
