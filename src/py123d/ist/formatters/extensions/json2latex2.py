#!/usr/bin/python
# -*- coding: utf-8 -*-
# author:   Jan Hybs
from py123d.ist.base import InputType
from py123d.ist.utils.texlist2 import TexList
from py123d.utils.logger import Logger
import json
from py123d.ist.extras import TypeReference

class LatexRecordDefault(object):
    """
    Class LatexRecordDefault is default formatter if no other formattor macthes
    """

    @staticmethod
    def default_format(tex, default):
        """
        :type tex: TexList
        :type default: ist.extras.TypeRecordKeyDefault
        """

        tex.append(tex._SPACE)
        tex._function_call(
            func=tex._IT,
            args=str(default.value).capitalize(),
            mode=tex.TYPE_PLAIN
        )
        #tex.append(tex._SPACE)

    @staticmethod
    def value_format(tex, default):
        """
        :type tex: TexList
        :type default: ist.extras.TypeRecordKeyDefault
        """

        tex.append(tex._SPACE)
        tex._function_call(
            func="ValueDefault",
            args=json.dumps(default.value),
            mode=tex.TYPE_PLAIN
        )
        #tex.append(tex._SPACE)

    @staticmethod
    def at_readtime_format(tex, default):
        """
        :type tex: TexList
        :type default: ist.extras.TypeRecordKeyDefault
        """
        tex.append('implicit value: "')
        tex.add(default.value, tex.TYPE_PLAIN)
        tex.append('"')

    @staticmethod
    def format(tex, default):
        """
        :type tex: TexList
        :type default: ist.extras.TypeRecordKeyDefault
        """
        LatexRecordDefault.format_rules.get(default.type).__func__(tex, default)

    format_rules = {
        'value at read time': at_readtime_format,
        'value at declaration': value_format,
        'optional': default_format,
        'obligatory': default_format,
        'default': default_format
    }


class LatexRecord(TexList):
    """
    Class LatexRecord is formatter class for type record
    """
    latex_name = 'RecordType'

    def format(self, record):
        """
        % begin{RecordType}
        %       {<href id>}
        %       {<record name>}                 % name of the record, used for header and for hypertarget in form IT::<record name>
        %       {<parent abstract record>}      % possible parent abstract record
        %       {<default conversion key>}      % possible auto conversion key
        %       {<link>}                        % possible hyperlink into hand written text
        %       {< record description>}         % description of the record
        %
        %       \KeyItem{<name>}                % name of the key
        %               {<type>}                % type of the key
        %               {<default value>}       % type of default value and possibly the value itself
        %               {<link>}                % possible hyperlink to hand written text
        %               {<key description>}     % description of the key
        %       ...
        % end{RecordType}

        :type record: ist.nodes.TypeRecord
        """
        self.begin(self.latex_name)

        # name
        self._newline()
        self._tab()
        #with self:
        self.add(record.href_id)
        self._newline()
        self._tab()
        self.add(record.href_name, self.TYPE_PLAIN)
            #self.macro_hyper_b(record)
        # implements
        self._newline()
        self._tab()
        with self:
            if record.implements:
                generic_roots=[abstract_ref.get_reference().get_generic_root().id for abstract_ref in record.implements]
                unique_roots=list(set(generic_roots))
                for impl in unique_roots:
                    self.macro_alink(TypeReference(impl).get_reference())
                    self.add(", ")
                self.pop()
        self.comment("implements")
        # conversion key
        self._newline()
        self._tab()
        with self:
            if record.reducible_to_key:
                self.macro_alink(record.reducible_to_key)
        self.comment("reducible to key")
        # hyperlink into hand written text TODO
        # LATER it can removed since is not used anymore
        #self._newline()
        #self._tab()
        #with self:
        #    pass
        #self.comment("OBSOLETE - hyperlink into hand written text")
        # description
        self._newline()
        self._tab()
        with self:
            self.append(self.description(record.description))
        # keys
        self._newline()
        for key in (record.keys or []):
            if key.href_name == "TYPE":
                continue
            self._newline()
            self._tab(2)
            with self.item_open('RecKey'):
                self.macro_key(key)

        self.end(self.latex_name)

    def macro_key(self, record_key):
        """
        :type record_key: ist.extras.TypeRecordKey
        """
        # name
        self._newline()
        self._tab(3)
        #with self:
        #   self.macro_hyper_b(record_key)
        self.add(record_key.href_id)
        self._newline()
        self._tab(3)
        self.add(record_key.href_name, self.TYPE_PLAIN)

        # type and parameters
        self._newline()
        self._tab(3)
        with self:
            ref = record_key.type.get_reference()
            parameters=self.get_key_type(ref)
        with self:
            self.format_parameters(parameters)

        # default
        self._newline()
        self._tab(3)
        with self:
            LatexRecordDefault.format(self, record_key.default)
        # hyperlink into hand written text TODO
        # LATER it can removed since is not used anymore
        #self._newline()
        #self._tab(3)
        #with self:
        #    pass
        #self.comment("OBSOLETE - hyperlink into hand written text")
        # description
        self._newline()
        self._tab(3)
        with self:
            self.append(self.description(record_key.description))

    def format_parameters(self, params):
        # params is SmartList
        if params:
            assert len(params) > 0
            for item in params:
                # item is TypeAttributeParameter
                param_type = item.reference.get_reference()
                self.add(item.name, self.TYPE_PLAIN)
                self.add(" = ")
                self.macro_alink(param_type)
                self.add(", ")
            self.pop()

    def get_key_type(self, ref):
        if ref.input_type == InputType.MAIN_TYPE:
            # selection, record, abstract, tuple
            if ref.has_generic_link():
                self.add("gen. " + str(ref.input_type).lower() + ": ")
                self.macro_alink(ref.get_generic_root())
                param_dict = ref.get_parameter_dict()
                return param_dict
            else:
                self.add(str(ref.input_type).lower() + ": ")
                # always point to generic root
                # if no generic type exists point to this reference
                self.macro_alink(ref.get_generic_root())
        elif ref.input_type == InputType.ARRAY:
            name = "array {range} of ".format(range=ref.range)
            self.add(name)
            # always point to generic root
            # if no generic type exists point to this reference
            return self.get_key_type(ref.subtype.target)
        elif ref.input_type == InputType.PARAMETER:
            name = "parameter: "
            self.add("parameter: "+ ref.name, self.TYPE_PLAIN)
        else:
            ref_range = (' ' + str(ref.get('range') or '')).rstrip()
            name = str(ref.input_type).capitalize()
            self.add(name + ref_range)
        return None

class LatexTuple(LatexRecord):
    """
    Class LatexTuple is formatter for type tuple
    """
    latex_name = 'TupleType'

    def __init__(self):
        super(LatexRecord, self).__init__()


class LatexSelection(TexList):
    """
    Class LatexSelection is formatter for type selection
    """

    def format(self, selection):
        """
        % begin{SelectionType}
        %       {<selection name>}
        %       {< selection description>}
        %
        %       \KeyItem
        %           {<value name>}
        %           {Key value description.}
        % end{SelectionType}

        :type selection: ist.nodes.TypeSelection
        """
        self.begin('SelectionType')

        # name
        self._newline()
        self._tab()
        #with self:
        #    self.macro_hyper_b(selection)
        self.add(selection.href_id)
        self._newline()
        self._tab()
        self.add(selection.href_name, self.TYPE_PLAIN)

        # description
        self._newline()
        self._tab()
        with self:
            self.append(self.description(selection.description))

        # values
        self._newline()
        for key in (selection.values or []):
            self._newline()
            self._tab(2)
            with self.item_open('SelectionItem'):
                self.macro_value(key)

        self.end('SelectionType')

    def macro_value(self, selection_value):
        """
        :type record_key: ist.extras.TypeSelectionValue
        """
        # name
        self._newline()
        self._tab(3)
        #with self:
        #    self.macro_hyper_b(selection_value)
        self.add(selection_value.href_id)
        self._newline()
        self._tab(3)
        self.add(selection_value.href_name, self.TYPE_PLAIN)
        # description
        self._newline()
        self._tab(3)
        with self:
            self.append(self.description(selection_value.description))


class LatexAbstractRecord(TexList):
    """
    Class LatexAbstractRecord is formatter for type abstract
    """

    def format(self, abstract_record):
        """
        % begin{AbstractType}
        %       {<record name>}
        %       {<default descendant>}
        %       {<link>}
        %       {<description>}
        %       \Descendant{<type name>}
        % end{AbstractType}
        :type abstract_record: ist.nodes.TypeAbstract
        """
        self.begin('AbstractType')

        # name
        self._newline()
        self._tab()
        #with self:
        #    self.macro_hyper_b(selection)
        self.add(abstract_record.href_id)
        self._newline()
        self._tab()
        self.add(abstract_record.href_name, self.TYPE_PLAIN)

        # descendant
        self._newline()
        self._tab()
        with self:
            if abstract_record.default_descendant:
                self.macro_alink(abstract_record.default_descendant.get_reference())
        #with self:
        #    self.macro_add_doc(abstract_record)
        # description
        self._newline()
        self._tab()
        with self:
            self.append(self.description(abstract_record.description))

        # implementations
        self._newline()
        for impl in abstract_record.implementations:
            self._newline()
            self._tab(2)
            with self.item_open('Descendant'):
                with self:
                    self.macro_alink(impl.get_reference())
            self.comment("") # end line comment to avoid spurious spaces
        self.end('AbstractType')


class LatexFormatter(object):
    """
    Class LatexFormatter is main formatter for entire document
    """

    formatters = {
        'TypeRecord': LatexRecord,
        # 'TypeRecordKey': LatexRecordKey,
        'TypeAbstract': LatexAbstractRecord,
        'TypeAbstractRecord': LatexAbstractRecord,
        # 'TypeString': LatexString,
        'TypeSelection': LatexSelection,
        'TypeTuple': LatexTuple,
        # 'TypeArray': LatexArray,
        # 'TypeInteger': LatexInteger,
        # 'TypeDouble': LatexDouble,
        # 'TypeBool': LatexBool,
        # 'TypeFilename': LatexFileName
    }

    @staticmethod
    def format(items):
        tex = TexList()
        #Logger.instance().info('Processing items...')
        for item in items:
            # do no format certain objects
            #Logger.instance().info('processing: %s (%s)' % (item, item.href_id))
            if not item.include_in_format():
                #Logger.instance().info('  - skipped')
                continue

            #Logger.instance().info('  +++ formatting +++')
            fmt = LatexFormatter.get_formatter_for(item)
            if fmt is not None:
                fmt.format(item)
                tex.extend(fmt)

        return tex

    @staticmethod
    def get_formatter_for(o):
        cls = LatexFormatter.formatters.get(o.__class__.__name__, None)
        if cls is None:
            return None
        return cls()
