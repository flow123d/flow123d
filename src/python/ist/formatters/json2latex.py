# encoding: utf-8
# author:   Jan Hybs


"""
\begin{RecordType}{\hyperB{IT::SequentialCoupling}{SequentialCoupling}}{\Alink{IT::Problem}{Problem}}{}{}{Record with data for a general sequential coupling.
}
\KeyItem{\hyperB{SequentialCoupling::TYPE}{TYPE}}{selection: Problem\_TYPE\_selection}{SequentialCoupling}{}{Sub-record selection.}
\KeyItem{\hyperB{SequentialCoupling::description}{description}}{String (generic)}{\textlangle{\it optional }\textrangle}{}{Short description of the solved problem.\\Is displayed in the main log, and possibly in other text output files.}
\KeyItem{\hyperB{SequentialCoupling::mesh}{mesh}}{record: \Alink{IT::Mesh}{Mesh}}{\textlangle{\it obligatory }\textrangle}{}{Computational mesh common to all equations.}
\KeyItem{\hyperB{SequentialCoupling::time}{time}}{record: \Alink{IT::TimeGovernor}{TimeGovernor}}{\textlangle{\it optional }\textrangle}{}{Simulation time frame and time step.}
\KeyItem{\hyperB{SequentialCoupling::primary-equation}{primary\_equation}}{abstract type: \Alink{IT::DarcyFlowMH}{DarcyFlowMH}}{\textlangle{\it obligatory }\textrangle}{}{Primary equation, have all data given.}
\KeyItem{\hyperB{SequentialCoupling::secondary-equation}{secondary\_equation}}{abstract type: \Alink{IT::Transport}{Transport}}{\textlangle{\it optional }\textrangle}{}{The equation that depends (the velocity field) on the result of the primary equation.}
\end{RecordType}
"""
from ist.nodes import Bool, AbstractRecord, Selection, Integer, DescriptionNode, String, ISTNode
from ist.utils.texlist import texlist

ur"""
\begin{RecordType}{\hyperB{IT::Root}{Root}}{}{}{}{Root record of JSON input for Flow123d.}
\KeyItem{\hyperB{Root::problem}{problem}}{abstract type: \Alink{IT::Problem}{Problem}}{\textlangle{\it obligatory }\textrangle}{}{Simulation problem to be solved.}
\KeyItem{\hyperB{Root::pause-after-run}{pause\_after\_run}}{Bool}{false}{}{If true, the program will wait for key press before it terminates.}
\end{RecordType}
"""


class LatexItemFormatter (object):
    def __init__ (self, tag_name=None):
        self.tag_name = tag_name

    def format (self, element):
        raise NotImplementedError ('Method format not implemented {}'.format (self.__class__.__name__))

    def format_as_child (self, *args):
        raise NotImplementedError ('Method format_as_child not implemented {}'.format (self.__class__.__name__))


class LatexSelection (LatexItemFormatter):
    def __init__ (self):
        super (LatexSelection, self).__init__ ('SelectionType')

    # "input_type": "Record",
    # "type_name": "Partition",
    # "type_full_name": "Partition",
    # "description": "Setting for various types of mesh partitioning.",
    # "reducible_to_key": "graph_type",
    # "keys": [
    # {
    # "key": "tool",
    # "type": "<<REFERENCE_ID>>",
    # "description": "Software package used for partitioning. See corresponding selection.",
    # "default": {
    # "type": "value at declaration",
    # "value": "METIS"
    # },
    # },
    #
    #
    # {selection: \Alink{IT::PartTool}{PartTool}}
    # {METIS}
    # {}
    # {Software package used for partitioning. See corresponding selection.}

    def format_as_child (self, self_selection, record_key, record):
        tex = texlist ()
        tex.KeyItem ()
        with tex:
            tex.hyperB (record_key.key, record.type_name)

        with tex:
            tex.append ('selection: ')
            if self_selection.include_in_format ():
                tex.Alink (self_selection.name)
            else:
                tex.append (tex.escape (self_selection.name))

        tex.add_s (record_key.default.value)
        tex.add ()
        tex.add_description_field (record_key.description)

        return tex

    # {
    # "id": "f9756fb2f66076a1",
    # "input_type": "Selection",
    # "name": "PartTool",
    # "full_name": "PartTool",
    # "description": "Select the partitioning tool to use.",
    # "values": [
    # {
    # "name": "PETSc",
    # "description": "Use PETSc interface to various partitioning tools."
    # },
    # {
    # "name": "METIS",
    # "description": "Use direct interface to Metis."
    # }
    # ]
    # }
    #
    # \begin{SelectionType}{\hyperB{IT::PartTool}{PartTool}}{Select the partitioning tool to use.}
    # \KeyItem{PETSc}{Use PETSc interface to various partitioning tools.}
    # \KeyItem{METIS}{Use direct interface to Metis.}
    # \end{SelectionType}
    #
    # SelectionType environment
    # usage:
    # \begin{SelectionType}{<selection name>}{< selection description>}
    # \KeyItem{<value name>}{<value>}
    # Key value description.
    # \end{SelectionType}

    def format (self, selection):
        tex = texlist (self.tag_name)

        with tex.element ():
            with tex:
                tex.hyperB (selection.name)
            tex.add_description_field (selection.description)

            for selection_value in selection.values:
                tex.newline ()
                tex.KeyItem (selection_value.name, selection_value.description)
            tex.newline ()

        return tex


# {
# "key": "tool",
# "description": "Software package used for partitioning. See corresponding selection.",
# "default": {
# "type": "value at declaration",
# "value": "METIS"
# },
# "type": "f9756fb2f66076a1"
# },
#
# \KeyItem{<name>}                % name of the key
# {<type>}                % type of the key
# {<default value>}       % type of default value and possibly the value itself
# {<link>}                %  possible hyperlink to hand written text
# {<key description>}     % description of the key
#
# \KeyItem{\hyperB{Partition::tool}{tool}}
# {selection: \Alink{IT::PartTool}{PartTool}}
# {METIS}
# {}
# {Software package used for partitioning. See corresponding selection.}


class LatexRecordKey (LatexItemFormatter):
    def __init__ (self):
        super (LatexRecordKey, self).__init__ ('KeyItem')


    def format (self, record_key, record):
        tex = texlist ()
        reference = record_key.type.get_reference ()

        # try to grab formatter and format type and default value based on reference type
        fmt = LatexFormatter.get_formatter_for (reference)
        if not fmt:
            print ' <<Missing formatter for {}>>'.format (type (reference))
        else:
            tex.extend (fmt.format_as_child (reference, record_key, record))

        return tex


# AbstractType environment
# usage:
# \begin{AbstractType}
# {<record name>}
# {<default descendant>}
# {<link>}
# {<description>}         % Description paragraph of the abstract type.
# \Descendant{<type name>}
# \end{AbstractType}

class LatexAbstractRecord (LatexItemFormatter):
    def __init__ (self):
        super (LatexAbstractRecord, self).__init__ ('AbstractType')

    # {
    # "key": "primary_equation",
    # "description": "Primary equation, have all data given.",
    # "default": {
    # "type": "obligatory",
    # "value": "OBLIGATORY"
    # },
    # "type": "89b3f40b6e805da8"
    # },
    # \KeyItem{\hyperB{SequentialCoupling::primary-equation}{primary\_equation}}
    # {abstract type: \Alink{IT::DarcyFlowMH}{DarcyFlowMH}}
    # {\textlangle{\it obligatory }\textrangle}
    # {}
    # {Primary equation, have all data given.}

    def format_as_child (self, abstract_record, record_key, record):
        tex = texlist ()
        tex.KeyItem ()
        with tex:
            tex.hyperB (record_key.key, record.type_name)

        with tex:
            tex.append ('abstract type: ')
            tex.Alink (abstract_record.name)

        tex.extend (
            LatexRecordKeyDefault ().format_as_child (record_key.default, record_key, record)
        )
        tex.add ()
        tex.add_description_field (record_key.description)

        return tex


    # {
    # "id": "89b3f40b6e805da8",
    # "input_type": "AbstractRecord",
    # "name": "DarcyFlowMH",
    # "full_name": "DarcyFlowMH",
    # "description": "Mixed-Hybrid  solver for saturated Darcy flow.",
    # "implementations": [
    # "3d0fe10b33e5936b",
    # "8d74b0bcfe5f8833",
    # "5bf2d1a105220256"
    # ]
    # },
    # \begin{AbstractType}
    # {\hyperB{IT::DarcyFlowMH}{DarcyFlowMH}}
    # {}
    # {}
    # {Mixed-Hybrid  solver for saturated Darcy flow.}
    # \Descendant{\Alink{IT::Steady-MH}{Steady\_MH}}
    # \Descendant{\Alink{IT::Unsteady-MH}{Unsteady\_MH}}
    # \Descendant{\Alink{IT::Unsteady-LMH}{Unsteady\_LMH}}
    # \end{AbstractType}


    def format (self, abstract_record):
        tex = texlist (self.tag_name)
        with tex.element ():
            with tex:
                tex.hyperB (abstract_record.name)

            if abstract_record.default_descendant:
                reference = abstract_record.default_descendant.get_reference ()
                with tex:
                    tex.Alink (reference.type_name)
                with tex:
                    tex.AddDoc (abstract_record.name)
            else:
                tex.add ()
                tex.add ()
            tex.add_description_field (abstract_record.description)

            for descendant in abstract_record.implementations:
                tex.newline ()
                tex.tag ('Descendant')
                with tex:
                    tex.Alink (descendant.get_reference ().type_name)
            tex.newline ()

        return tex

# {
# "id": "29b5533100b6f60f",
# "input_type": "String",
# "name": "String",
# "full_name": "String"
# }
#
# %       \KeyItem{<name>}                % name of the key
# %               {<type>}                % type of the key
# %               {<default value>}       % type of default value and possibly the value itself
# %               {<link>}                % possible hyperlink to hand written text
# %               {<key description>}     % description of the key
#
# {String (generic)}
# {\textlangle{\it optional }\textrangle}
# {}
# {description}


class LatexString (LatexItemFormatter):
    def __init__ (self):
        super (LatexString, self).__init__ ('String')


    # {
    # "id": "b9614d55a6c3462e",
    # "input_type": "Record",
    # "type_name": "Region",
    # "type_full_name": "Region",
    # "description": "Definition of region of elements.",
    # "keys": [
    # {
    # "key": "name",
    # "type": "<<REFERENCE_ID>>",
    # "description": "Label (name) of the region. Has to be unique in one mesh.\n",
    # "default": {
    # "type": "obligatory",
    # "value": "OBLIGATORY"
    # },
    #
    # \KeyItem{\hyperB{Region::name}{name}}
    # {String (generic)}
    # {\textlangle{\it obligatory }\textrangle}
    # {}
    # {Label (name) of the region. Has to be unique in one mesh.\\}

    def format_as_child (self, self_string, record_key, record):
        tex = texlist ()
        tex.KeyItem ()
        with tex:
            tex.hyperB (record_key.key, record.type_name)
        tex.add ('String (generic)')

        with tex:
            tex.textlangle (record_key.default.type)
        tex.add ()
        tex.add_description_field (record_key.description)

        return tex


# RecordType environment
#
# usage:
# \begin{RecordType}
# {<record name>}                 % name of the record, used for header and for hypertarget in form IT::<record name>
# {<parent abstract record>}      % possible parent abstract record
# {<default conversion key>}      % possible auto conversion key
# {<link>}                        % possible hyperlink into hand written text
# {< record description>}         % description of the record
#
# \KeyItem{<name>}                % name of the key
# {<type>}                % type of the key
# {<default value>}       % type of default value and possibly the value itself
# {<link>}                %  possible hyperlink to hand written text
# {<key description>}     % description of the key
# ...
# \end{RecordType}

class LatexRecord (LatexItemFormatter):
    def __init__ (self):
        super (LatexRecord, self).__init__ ('RecordType')


    # {
    # "key": "mesh",
    # "description": "Computational mesh common to all equations.",
    # "default": {
    # "type": "obligatory",
    # "value": "OBLIGATORY"
    # },
    # "type": "c57e1ac33a446313"
    # },
    #
    # \KeyItem{\hyperB{SequentialCoupling::mesh}{mesh}}
    # {record: \Alink{IT::Mesh}{Mesh}}
    # {\textlangle{\it obligatory }\textrangle}
    # {}
    # {Computational mesh common to all equations.}


    def format_as_child (self, self_record, record_key, record):
        tex = texlist ()
        tex.KeyItem ()
        with tex:
            tex.hyperB (record_key.key, record.type_name)

        with tex:
            tex.append ('record: ')
            tex.Alink (self_record.type_name)

        tex.extend (
            LatexRecordKeyDefault ().format_as_child (record_key.default, record_key, record)
        )
        tex.add ()
        tex.add_description_field (record_key.description)

        return tex

    # {
    # "id": "81a9cc0d6917a75",
    # "input_type": "Record",
    # "type_name": "SequentialCoupling",
    # "type_full_name": "SequentialCoupling:Problem",
    # "description": "Record with data for a general sequential coupling.\n",
    # "implements": [
    #      "1b711c7cb758740"
    #   ],
    #   "keys": [
    #      {
    #         "key": "TYPE",
    #         "description": "Sub-record selection.",
    #         "default": {
    #            "type": "value at declaration",
    #            "value": "SequentialCoupling"
    #         },
    #         "type": "b0bf265898e2625b"
    #      },
    #      ...
    #      {
    #         "key": "secondary_equation",
    #         "description": "The equation that depends (the velocity field) on the result of the primary equation.",
    #         "default": {
    #            "type": "optional",
    #            "value": "OPTIONAL"
    #         },
    #         "type": "ba303ae783dcf903"
    #      }
    #   ]
    # }
    # \begin{RecordType}
    #       {\hyperB{IT::SequentialCoupling}{SequentialCoupling}}
    #       {\Alink{IT::Problem}{Problem}}
    #       {}
    #       {}
    #       {Record with data for a general sequential coupling.}
    # \KeyItem{\hyperB{SequentialCoupling::TYPE}{TYPE}}
    #     {selection: Problem\_TYPE\_selection}{SequentialCoupling}
    #     {}
    #     {Sub-record selection.}
    # ...
    # \KeyItem{\hyperB{SequentialCoupling::secondary-equation}{secondary\_equation}}
    #     {abstract type: \Alink{IT::Transport}{Transport}}
    #     {\textlangle{\it optional }\textrangle}
    #     {}
    #     {The equation that depends (the velocity field) on the result of the primary equation.}
    # \end{RecordType}


    def format (self, record):
        tex = texlist (self.tag_name)
        reference_list = record.implements


        print record, type(record)
        try: print record.keys
        except Exception as e: print e
    
        with tex.element ():
            with tex:
                tex.hyperB (record.type_name)

            # TODO what if multiple inheritance? list
            if reference_list:
                with tex:
                    for reference in reference_list:
                        tex.Alink (reference.get_reference ().name)
            else:
                tex.add ()

            if record.reducible_to_key:
                with tex:
                    tex.Alink (record.reducible_to_key, record.type_name)
            else:
                tex.add ()

            # hyperlink into hand written text
            # LATER it can removed since is not used anymore
            tex.add ()
            tex.add_description_field (record.description)

            # record keys
            # for record_key in record.keys:
            #     tex.newline ()
            #     fmt = LatexFormatter.get_formatter_for (record_key)
            #     tex.extend (fmt.format (record_key, record))
            tex.newline ()

        return tex


class LatexRecordKeyDefault (LatexItemFormatter):
    def __init__ (self):
        super (LatexRecordKeyDefault, self).__init__ ('')

        self.format_rules = {
            'value at read time': self.raw_format,
            'value at declaration': self.textlangle_format,
            'optional': self.textlangle_format,
            'obligatory': self.textlangle_format
        }

    def format_as_child (self, self_default, record_key, record):
        method = self.format_rules.get (self_default.type, None)
        if method:
            return method (self_default, record_key, record)

        return LatexRecordKeyDefault.textlangle_format (self_default, record_key, record)

    def textlangle_format (self, self_default, record_key, record):
        tex = texlist ()
        with tex:
            tex.textlangle (self_default.value)
        return tex

    def raw_format (self, self_default, record_key, record):
        tex = texlist ()
        tex.add_s (self_default.value)
        return tex


class LatexUniversal (LatexItemFormatter):
    def __init__ (self):
        super (LatexUniversal, self).__init__ ('')

    # \KeyItem
    # {\hyperB{Mesh::regions}{regions}}
    # {Array  of record: \Alink{IT::Region}{Region}}
    # {\textlangle{\it optional }\textrangle}
    # {}
    # {List of additional region definitions not contained in the mesh.}
    def _start_format_as_child (self, self_object, record_key, record):
        tex = texlist ()
        tex.KeyItem ()
        with tex:
            tex.hyperB (record_key.key, record.type_name)
        return tex

    def _end_format_as_child (self, self_object, record_key, record):
        tex = texlist ()

        tex.extend (
            LatexRecordKeyDefault ().format_as_child (record_key.default, record_key, record)
        )
        tex.add ()
        tex.add_description_field (record_key.description)

        return tex

    def _format_as_child (self, self_object, record_key, record):
        raise Exception ('Not implemented yet')

    def format_as_child (self, self_object, record_key, record):
        tex = texlist ()
        tex.extend (self._start_format_as_child (self_object, record_key, record))
        tex.extend (self._format_as_child (self_object, record_key, record))
        tex.extend (self._end_format_as_child (self_object, record_key, record))
        return tex


class LatexArray (LatexUniversal):
    def _format_as_child (self, self_array, record_key, record):
        subtype = self_array.subtype.get_reference ()
        tex = texlist ()
        with tex:
            if type (subtype) == Integer:
                tex.append ('Array of {subtype} {subrange}'.format (
                    range=self_array.range, subtype=subtype.input_type,
                    subrange=subtype.range))
            else:
                tex.append ('Array{range} of {subtype}'.format (
                    range=' ' + str (self_array.range) if not self_array.range.is_pointless () else '',
                    subtype=subtype.input_type))

            if type (subtype) == String:
                tex.append (' (generic)')

            if issubclass (subtype.__class__, DescriptionNode):
                tex.append (': ')
                tex.Alink (subtype.get ('type_name', 'name'))
            else:
                # no link
                pass

        return tex


class LatexInteger (LatexUniversal):
    def _format_as_child (self, self_int, record_key, record):
        tex = texlist ()
        tex.add ('Integer' + str (self_int.range))
        return tex

    def _end_format_as_child (self, self_object, record_key, record):
        tex = texlist ()
        tex.extend (
            LatexRecordKeyDefault ().format_as_child (record_key.default, record_key, record)
        )
        tex.add ()
        tex.add_description_field (record_key.description)
        return tex


class LatexDouble (LatexUniversal):
    def _format_as_child (self, self_double, record_key, record):
        tex = texlist ()
        tex.add ('Double' + str (self_double.range))
        return tex

    def _end_format_as_child (self, self_object, record_key, record):
        tex = texlist ()
        tex.extend (
            LatexRecordKeyDefault ().format_as_child (record_key.default, record_key, record)
        )
        tex.add ()
        tex.add_description_field (record_key.description)
        return tex


class LatexBool (LatexUniversal):
    def _format_as_child (self, self_bool, record_key, record):
        tex = texlist ()
        tex.add ('Bool')
        return tex

    def _end_format_as_child (self, self_object, record_key, record):
        tex = texlist ()
        tex.add (record_key.default.value) # todo LatexRecordKeyDefault
        tex.add ()
        tex.add_description_field (record_key.description)

        return tex


class LatexFileName (LatexUniversal):
    def _format_as_child (self, self_fn, record_key, record):
        tex = texlist ()
        tex.add (self_fn.file_mode + ' file name')
        return tex


class LatexFormatter (object):
    formatters = {
        'Record': LatexRecord,
        'RecordKey': LatexRecordKey,
        'AbstractRecord': LatexAbstractRecord,
        'String': LatexString,
        'Selection': LatexSelection,
        'Array': LatexArray,
        'Integer': LatexInteger,
        'Double': LatexDouble,
        'Bool': LatexBool,
        'FileName': LatexFileName
    }

    @staticmethod
    def format (items):
        tex = texlist ()

        for item in items:
            # format only IST nodes
            if issubclass (item.__class__, ISTNode):

                # do no format certain objects
                if not item.include_in_format ():
                    continue

                try:
                    fmt = LatexFormatter.get_formatter_for (item)
                    tex.extend (fmt.format (item))
                    tex.newline ()
                    tex.newline ()
                except NotImplementedError as e:
                    print e

        return tex

    @staticmethod
    def get_formatter_for (o):
        cls = LatexFormatter.formatters.get (o.__class__.__name__, None)
        if cls is None:
            return None
        return cls ()


        # TODO default_descendant