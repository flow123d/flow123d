'''
New convertor of input YAML files.
Features:
- conversion rules implemented directly in Python using predefined set of actions in form of methods
- conversion rules operates on YAML tree using natural dict and list composed types
- use ruamel.yaml lib to preserve comments, order of keys etc.
- each set of rules will be in separate method, registered into main list of rules,
  every such change set will have an unique number (increasing), some change sets may be noted by flow123d release
  Both the release and input change set number will be part of the YAML file in order to apply changes only once.
- Try to make actions reversible, so we can make also (some) back conversion.  
- can run in quiet mode or in debug mode, individual change sets may be noted as stable to do not report warnings as default
- can be applied to the input format specification and check that it produce target format specification

TODO:
- force the ruaml.yaml to read custom tags for all types, without complaining (add_mluit_constructor works)
- force to keep these tags in created objects

'''
import ruamel.yaml as ruml
import re
import os
#import sys
import argparse
import fnmatch
import types
import logging
#import copy


MapType=ruml.comments.CommentedMap
ListType=list


CommentsTag = ruml.comments.Tag

class CommentedScalar:
    """
    Class to store all scalars with their tags
    """
    original_constructors={}


    def __repr__(self):
        str(self.value)

    @classmethod
    def to_yaml(cls, dumper, data):
        representer = dumper.yaml_representers[type(data.value).__mro__[0]]
        node = representer(dumper, data.value)
        if data.tag.value is None:
            tag = node.tag
        elif data.tag.value.startswith(u'tag:yaml.org,2002'):
            tag = node.tag
        else:
            tag = data.tag.value
        print("val: ", data.value, "repr: ", node.value, "tag: ", tag)
        return dumper.represent_scalar(tag, node.value)


    def __init__(self, tag, value):
        self.tag.value = tag
        self.value = value

    @property
    def tag(self):
        # type: () -> Any
        if not hasattr(self, CommentsTag.attrib):
            setattr(self, CommentsTag.attrib, CommentsTag())
        return getattr(self, CommentsTag.attrib)



def construct_any_tag(self, tag_suffix, node):
    if tag_suffix is None:
        orig_tag = None
    else:
        orig_tag = "!" + tag_suffix
    if isinstance(node, ruml.ScalarNode):

        implicit_tag = self.composer.resolver.resolve(ruml.ScalarNode, node.value, (True, None))
        if implicit_tag in self.yaml_constructors:
            #constructor = CommentedScalar.original_constructors[implicit_tag]
            constructor = self.yaml_constructors[implicit_tag]
        else:
            constructor = self.construct_undefined

        data = constructor(self, node)
        if isinstance(data, types.GeneratorType):
            generator = data
            data = next(generator)  # type: ignore

        scal = CommentedScalar(orig_tag, data)
        yield scal

    elif isinstance(node, ruml.SequenceNode):
        for seq in self.construct_yaml_seq(node):
            seq.yaml_set_tag(orig_tag)
            yield seq
    elif isinstance(node, ruml.MappingNode):
        for map in self.construct_yaml_map(node):
            map.yaml_set_tag(orig_tag)
            yield map
    else:
        for dummy in self.construct_undefined(node):
            yield dummy
"""
def construct_scalar(self, node):
    gen = construct_any_tag(self, None, node)
    for item in gen:
        yield item
"""

#def dump_commented_scalar(cls, data):
#    data.dump(cls)

def represent_commented_seq(cls, data):
    if data.tag.value is None:
        tag = u'tag:yaml.org,2002:seq'
    else:
        tag = data.tag.value
    return cls.represent_sequence(tag, data)


yml=ruml.YAML(typ='rt')
yml.representer.add_representer(CommentedScalar, CommentedScalar.to_yaml)
yml.representer.add_representer(ruml.comments.CommentedSeq, represent_commented_seq)
yml.constructor.add_multi_constructor("!", construct_any_tag)
"""
CommentedScalar.original_constructors = copy.deepcopy(yml.constructor.yaml_constructors)

yml.constructor.add_constructor(
    u'tag:yaml.org,2002:null',
    construct_scalar)

yml.constructor.add_constructor(
    u'tag:yaml.org,2002:bool',
    construct_scalar)

yml.constructor.add_constructor(
    u'tag:yaml.org,2002:int',
    construct_scalar)

yml.constructor.add_constructor(
    u'tag:yaml.org,2002:float',
    construct_scalar)

yml.constructor.add_constructor(
    u'tag:yaml.org,2002:binary',
    construct_scalar)

yml.constructor.add_constructor(
    u'tag:yaml.org,2002:timestamp',
    construct_scalar)
yml.constructor.add_constructor(
    u'tag:yaml.org,2002:str',
    construct_scalar)

"""


'''
Helpers:
'''
def parse_yaml_str(yaml_string):
    '''
    Parse given string in YAML format and return the YAML tree.
    In combination with other methodsAdd key 'key_name' to the map at 'path', and assign to it the given 'key_value'.
    The 'key_value' can be scalar, dict or list.
    '''
    pass

def is_list_node(node):
    return type(node) in [ list, ruml.comments.CommentedSeq ]

def is_map_node(node):
    return type(node) in [ dict, ruml.comments.CommentedMap ]

def is_scalar_node(node):
    return type(node) in [ CommentedScalar, int, float, bool, None, str ]


class Changes:    
    def __init__(self):
        self.current_version=None
        self._version_changes=None
        self._changes=[]

    def new_version(self, version):
        """
        Close list of changes form previous version, open new change list.
        :param version: Version tag, e.g. "major.minor. patch"
        :return: None
        """
        if self.current_version :
            self._changes.append( (self.current_version, self._version_changes) )
        self.current_version=version
        self._version_changes=[]

    def __iadd__(self, change):
        """
        :param change:
        :return:
        """
        self._version_changes.append(change)
        return self



    def unify_tree_scalars(self, node):
        if is_list_node(node):
            for idx in range(len(node)):
                node[idx] = self.unify_tree_scalars(node[idx])
        elif is_map_node(node):
            for key, child in node.items():
                node[key] = self.unify_tree_scalars(child)
        elif is_scalar_node(node):
            if type(node) != CommentedScalar:

                tags_for_types = { \
                    float : u'tag:yaml.org,2002:float',\
                    int : u'tag:yaml.org,2002:int',\
                    bool : u'tag:yaml.org,2002:bool',\
                    str : u'tag:yaml.org,2002:str',\
                    None : u'tag:yaml.org,2002:null'}
                assert type(node) in tags_for_types
                tag = tags_for_types[type(node)]
                node = CommentedScalar(tag, node)
        else:
            assert False, "Unsupported node type: {}".format(type(node))
        return node

    def unify_trees(self, trees):
        return [ (fname, self.unify_tree_scalars(tree)) for fname, tree in trees ]

    def apply_changes(self, trees, in_version, out_version,  reversed=False, warn=True):
        """
        Apply initailized list of actions to the data tree 'root'.
        :param root: Input data tree.
        :param warn: Produce wranings for nondeterministic conversions.
        :param reversed: Produce backward conversion from the target version to the initial version.
        :return: Data tree after conversion.
        """
        trees = self.unify_trees(trees)

        # close last list
        self.new_version(None)

        # reverse
        if reversed:
            self._changes.reverse()
            for version, list in self._changes:
                list.reverse()

        active = False
        for version, change_list in self._changes:
            if version == in_version:
                active = True
            if version == out_version:
                break
            if active:
                for path_set, forward, backward, ac_name in change_list:
                    if reversed:
                        action=backward
                    else:
                        action=forward
                    changed_files=[]
                    for fname, tree in trees:
                        is_changed = path_set.apply(action, tree)
                        if is_changed:
                            changed_files.append(fname)
                    if changed_files:
                        print("Action: ", ac_name)
                        print("Applied to: ", changed_files)




class PathSet(object):
    """
    Set of places in YAML file to which apply a given rule
    """
    def __init__(self, path_patterns, **kwds):
        """
        Initialize the set by all plases of the given root structure using deep first search.
        :param path: Path to root (e.g. "/key1/0/key2/1"
        :param root: Root map or array.
        """
        self.patterns=[]
        for p in path_patterns:
            " ** = any number of levels, any key or index per level"
            pp = re.sub('\*\*', '[a-zA-Z0-9_]@(/[a-zA-Z0-9_]@)@',p)
            " * = single level, any key or index per level"
            pp = re.sub('\*', '[a-zA-Z0-9_]@', pp)
            " # = single level, only indices "
            pp = re.sub('\#', '[0-9]@', pp)
            pp = re.sub('@', '*', pp)
            pp = "^" + pp + "$"
            self.patterns.append(pp)
        self.options=kwds

    def apply(self, action, root):
        """
        Apply 'fun' to all paths in path set in data tree 'root'.
        :param fun: Function accepting a data tree path as parameter.
        :param root: Data tree
        :return: None
        """
        logging.debug("Patterns: {}".format(self.patterns))
        return self.dfs_apply(action, [root], "/")

    def dfs_apply(self, action, data_path, path):
        matches=[]
        current = data_path[-1]
        logging.debug("DFS apply: " + str(path))
        if self.match(data_path, path):
            matches+=[current]
            self.current_path=path
            action(data_path)
        if is_list_node(current):
            for idx, child in enumerate(current):
                matches += self.dfs_apply(action, data_path + [ child ], path  + str(idx) + "/"  )
        elif is_map_node(current):
            for key, child in current.items():
                matches += self.dfs_apply(action, data_path + [ child ] , path  + str(key) + "/" )
        return matches


    def match(self, data_path, path):
        for pattern in self.patterns:
            if re.match(pattern, path) != None:
                logging.debug("Match path")
                for key, param in self.options.items():
                    logging.debug("Match path")
                    if key == "have_tag":
                        tag_path = list(param.split('/'))
                        target = self.traverse_tree(data_path, tag_path[0:-1])
                        target = target[-1]
                        if not target or not hasattr(target, "tag"):
                            return False
                        if target.tag.value != "!" + tag_path[-1]:
                            return False
                logging.debug("Full Match")
                return True
        return False


    def traverse_tree(self, data_path, rel_path):
        """
        Move accross data tree starting at address given by 'data_path', according
        to rel_path.
        :param self:
        :param data_path: List of data nodes starting from root down to current node.
        :param rel_path: Relative address of the target node, e.g. "../key_name/1"
        :return: Data path of target node. None in case of incomaptible address.
        """
        target_path=data_path.copy()
        for key in rel_path:
            curr = target_path[-1]
            if key == "..":
                target_path.pop()
            elif key.isdigit():
                assert(type(curr) == ListType)
                idx = int(key)
                if idx >= len(curr):
                    return None
                target_path.append( curr[idx] )
            else:
                assert (type(curr) == MapType)
                if not key in curr:
                    return None
                target_path.append(curr[key])
        return target_path

    def __getattribute__(self, item):
        """
        For every action method (with underscore) is defined a proxy method to add the action into
        change list. Both forward and backward actions are generated.
        :param item:
        :return:
        """
        func_name="_"+item
        try:
            func = object.__getattribute__(self, func_name)
        except AttributeError:
            func=None
        if func:
            def wrap(*args, **kargs):
                def forward(path):
                    func(path, *args, reversed=False, **kargs)
                def backward(path):
                    func(path, *args, reversed=True, **kargs)
                return (self, forward, backward, func.__name__)
            return wrap
        else:
            return object.__getattribute__(self, item)


    '''
    Actions:
    Need a decorator that for the particular function creates a method to save the action.
    '''
    def _add_key_to_map(self, path, key, value, reversed):
        '''
        ACTION.
        For every path P in path set add key 'key_name' to the map at path P.
        This path must be a map. Assign the given 'key_value' to the key.
        The 'key_value' can be scalar, dict or list.
        REVERSE.
        For every path P in the path set 'path' remove key 'key_name' from the map.
        '''
        curr=path[-1]
        if not is_map_node(curr):
            logging.warning("Expecting map at path: {}".format(self.current_path))
            return
        if reversed:
            if key in curr:
                del curr[key]
        else:
            assert(not key in curr)
            for idx, dict_key in enumerate(curr.keys()):
                if str(dict_key) > str(key):
                    break
            curr.insert(idx, key, value)


    """
    def remove_key(path_set, key_name, key_value):
        '''
        ACTION.
        For every path P in path list 'path_set' remove key 'key_name'. This is an inverse action to add_key_to_map.
        REVERSE.
        For every path P in the path set 'path' add key 'key_name' to the map and set it to the given default value
        `key_value`.
        '''
    """


    def _rename_key(self, path, old_key, new_key, reversed):
        '''
        ACTION.
        For every path P in the path list 'path_set', which has to be a map.
        Rename its key 'old_name' to 'new_name'.
        REVERSE.
        For every path P in the path set 'path', rename 'new_name' to 'old_name'.
        TODO: generalize to move_key, problem where to apply reversed action.
        '''
        curr=path[-1]
        if not is_map_node(curr):
            logging.warning("Expecting map at path: {}".format(self.current_path))
            return
        if reversed:
            curr[old_key] = curr.pop(new_key)
        else:
            curr[new_key] = curr.pop(old_key)

    def _scale_scalar(self,  path, multiplicator, reversed):
        '''
        ACTION.
        For every path P in the path set 'path' which has to be a scalar, multiply it by 'multiplicator'
        REVERSE.
        For every path P in the path set 'path' which has to be a scalar, divide it by 'multiplicator'
        '''
        curr=path[-1]
        if not is_scalar_node(curr):
            #print("Node: ", curr)
            logging.warning("Expecting scalar at path: {} get: {}".format(self.current_path, type(curr)))
            return
        assert type(curr) == CommentedScalar, "Wrong type: {} {}".format(type(curr))
        if not type(curr.value) in [int, float]:
            return
        is_int = curr.value is int
        if reversed:
            curr.value /= multiplicator
        else:
            curr.value *= multiplicator
        if is_int and float(curr.value).is_integer():
            curr.value = int(curr.value)

    def _replace_value(self, path, re_forward, re_backward):
        """
        ACTION.
        For 'path' which has to be string scalar, apply regexp substitution.
        :param path:
        :param re_forward:  (regexp, substitute)
        :param re_backward: (regexp, substitute)
            ... used as re.sub(regexp, substitute, value)

        If regexp is None, then substitute is tuple for the manual conversion.
        :return: None
        """
        curr=path[-1]
        if not is_scalar_node(curr):
            #print("Node: ", curr)
            logging.warning("Expecting scalar at path: {} get: {}".format(self.current_path, type(curr)))
            return
        assert type(curr) == CommentedScalar, "Wrong type: {} {}".format(type(curr))
        if not curr.value is str:
            #print("Node: ", curr)
            logging.warning("Expecting string at path: {} get: {}".format(self.current_path, type(curr)))
            return
        if reversed:
            if re_backward[0] is None:
                pass    # manual
            else:
                re.sub(re_backward[0], re_backward[1], curr.value)
        else:
            if re_forward[0] is None:
                pass    # manual
            else:
                re.sub(re_forward[0], re_forward[1], curr.value)

    """
    def manual_conversion(self,  invalidate="key", message_forward, message_backward):
        '''
        ACTION.
        For every path P in the path set 'path' which has to end by key. Rename the key (if invalidate='key')
        of the tag (if invalidate='tag') by postfix '_NEED_EDIT'. And appended comment with the message_forward.
        REVERSE.
        For every path P in the path set 'path', make the same, but use message_backward for the comment.
        '''
        pass
    """

"""
Implement 'has_key', make tests of this class.
"""

changes = Changes()    
    
    
    
# Add header key 'flow123d_version'
changes.new_version("1.8.2")
changes += PathSet(["/"]).add_key_to_map(key = "flow123d_version", value = "2.0.0")

# Change degree keys in PadeApproximant
path_set = PathSet(["/problem/secondary_equation/**/ode_solver/"],
                   have_tag="PadeApproximant")
changes += path_set.rename_key(old_key="denominator_degree", new_key="pade_denominator_degree")
changes += path_set.rename_key(old_key="nominator_degree", new_key="pade_nominator_degree")

# Change sign of boundary fluxes
path_set = PathSet([
    "/problem/secondary_equation/input_fields/*/bc_flux/",
    "/problem/primary_equation/input_fields/*/bc_flux/",
    "/problem/secondary_equation/input_fields/*/bc_flux/#/",
    "/problem/primary_equation/input_fields/*/bc_flux/#/"])
changes += path_set.scale_scalar(multiplicator=-1)

path_set = PathSet([
    "/problem/secondary_equation/input_fields/*/bc_flux/value/",
    "/problem/primary_equation/input_fields/*/bc_flux/value/"],
    have_tag="../FieldConstant")
changes += path_set.scale_scalar(multiplicator=-1)

"""
path_set = PathSet([ "/problem/secondary_equation/input_fields/*/bc_flux/value/",
             "/problem/primary_equation/input_fields/*/bc_flux/value/"],
             have_tag="../FieldFormula")

changes += path_set.replace_value(
            re_forward=("^(.*)$", "-(\\1)"),
            re_backward=("^(.*)$", "-(\\1)"))

path_set = PathSet(["/problem/secondary_equation/input_fields/*/bc_flux",
             "/problem/primary_equation/input_fields/*/bc_flux"],
                   have_tag="../FieldElementwise")

changes += path_set.manual_conversion(invalidate='tag',
                             message_forward="Change sign of this field in the GMSH file.",
                             message_backward="Change sign of this field in the GMSH file.")

# Replace Robin and Neumann conditions by total flux
path = PathSet(["/problem/primary_equation/input_fields/*/bc_type"])
changes+= path.replace_value(path, old_value=["neumann", "robin"], new_value="total_flux", manual_reversion="manual")

path = [ PathPattern("/problem/secondary_equation/input_fields/*/bc_type"]
changes+= change_value(path, old_value=["neumann", "robin"], new_value="diffusive_flux", manual_reversion="manual")  
                         

# Rename equations and couplings
changes.new_set()

path = "/problem/secondary_equation"
changes += rename_tag(path, old_tag="TransportOperatorSplitting", new_tag="Coupling_OperatorSplitting")

path_os = [path, HaveTag("Coupling_OperatorSplitting")]
destination_path="/problem/secondary_equation/transport"
changes += add_key( path_os, key="transport", tag="Solute_Advection_FV") # parent="Transport" #seems unnecessary
changes += move_keys( path_os, keys=["output_fields", "input_fields"], destination_path )

path_dg = [path, HaveTag("SoluteTransport_DG")]
changes += add_key( path_dg, key="transport", tag="Solute_AdvectionDiffusion_DG")
changes += move_keys( path_dg, keys=["input_fields", "output_fields", "solver", "dg_order", "dg_variant", "solvent_density"], destination_path)
changes += rename_tag(path, old_tag="SoluteTransport_DG", new_tag="Coupling_OperatorSplitting")

changes += rename_tag(path,  old_tag="HeatTransfer_DG", new_tag="Heat_AdvectionDiffusion_DG")

# Remove r_set, use region instead
changes.new_set()
changes += move_and_rename_key(PathPattern("/**/input_fields/*/r_set"), new_path="$1/input_fields/$2/region")

# Changes in mesh record
changes.new_set()
path = PathPattern("/problem/mesh/regions/*/")
destination_path = "/problem/mesh/_regions_elementary/$1/"
changes += move_keys(path, keys=["name","id", "element_list"], destination_path)
path_el = [ PathPattern("/problem/mesh/_regions_elementary/$1/"), HaveKey("element_list") ]
changes += set_tag(path_el, new_tag="From_Elements")
path_el = [ PathPattern("/problem/mesh/_regions_elementary/$1/"), HaveKey("id") ]
changes += set_tag(path_el, new_tag="From_ID")

path = PathPattern("/problem/mesh/sets/*/")
destination_path = "/problem/mesh/_sets_new/$1/"
changes += move_keys(path, keys=["name","region_ids", "region_labels", ], destination_path)
path_el = [ PathPattern("/problem/mesh/_regions_elementary/$1/"), HaveKey("element_list") ]
changes += set_tag(path_el, new_tag="From_Elements")
path_el = [ PathPattern("/problem/mesh/_regions_elementary/$1/"), HaveKey("id") ]
changes += set_tag(path_el, new_tag="From_ID")
"""

'''
     { 
      "NAME" : "mesh sets setup, name",
      "action": "move-key",
      "parameters": {
        "source_path":"/problem/mesh/sets/*/name",
        "destination_path":"/problem/mesh/sets_new/$1/name",
          "create_path":true
      }
    },   
    { 
      "NAME" : "mesh sets setup, elementary regions setup, region_ids",
      "action": "move-key",
      "parameters": {
        "source_path":"/problem/mesh/sets/*/region_ids",
        "destination_path":"/problem/mesh/sets_new/$1/region_ids",
        "set_type_path":"/problem/mesh/sets_new/$1",
        "new_type":"Union",
        "create_path":true  
      }
    },
    { 
      "NAME" : "mesh sets setup, elementary regions setup, region_labels",
      "action": "move-key",
      "parameters": {
        "source_path":"/problem/mesh/sets/*/region_labels",
        "destination_path":"/problem/mesh/sets_new/$1/regions",
        "set_type_path":"/problem/mesh/sets_new/$1",
        "new_type":"Union", 
        "create_path":true 
      }
    },
    { 
      "NAME" : "mesh sets setup, elementary regions setup, union",
      "action": "move-key",
      "parameters": {
        "source_path":"/problem/mesh/sets/*/union",
        "destination_path":"/problem/mesh/sets_new/$1/regions",
        "set_type_path":"/problem/mesh/sets_new/$1",
        "new_type":"Union",
        "create_path":true  
      }
    },
    { 
      "NAME" : "mesh sets setup, elementary regions setup, intersection",
      "action": "move-key",
      "parameters": {
        "source_path":"/problem/mesh/sets/*/intersection",
        "destination_path":"/problem/mesh/sets_new/$1/regions",
        "set_type_path":"/problem/mesh/sets_new/$1",
        "new_type":"Intersection",
        "create_path":true  
      }
    },
    { 
      "NAME" : "mesh sets setup, elementary regions setup, difference",
      "action": "move-key",
      "parameters": {
        "source_path":"/problem/mesh/sets/*/difference",
        "destination_path":"/problem/mesh/sets_new/$1/regions",
        "set_type_path":"/problem/mesh/sets_new/$1",
        "new_type":"Difference",
        "create_path":true  
      }
    },
    {
      "action": "delete-key",
      "parameters": {
        "path": "/problem/mesh/sets/*",
        "deep": false   
      }
    },
    {
      "action": "delete-key",
      "parameters": {
        "path": "/problem/mesh/sets",
        "deep": false   
      }
    },
    {
      "action": "delete-key",
      "parameters": {
        "path": "/problem/mesh/regions/*",
        "deep": false
      }
    },
    {
      "action": "delete-key",
      "parameters": {
        "path": "/problem/mesh/regions",
        "deep": false
      }
    },
    {
      "action" : "merge-arrays",
      "parameters": {
        "source_path":"/problem/mesh/regions_elementary",
        "addition_path":"/problem/mesh/sets_new",
        "destination_path":"/problem/mesh/regions"        
      }
    },           
    {
      "action" : "merge-arrays",
      "parameters": {
        "source_path":"/problem/mesh/sets_new",
        "addition_path":"/problem/mesh/regions_elementary",
        "destination_path":"/problem/mesh/regions"        
      }
    },


    {
      "NAME": "Move linear_solver under nonlinear_solver in DarcyFlow.",
      "action": "add-key",
      "parameters": {
        "path": "/problem/primary_equation",
        "key": "nonlinear_solver"
      }
    }, 
    {
      "NAME": "Move linear_solver under nonlinear_solver in DarcyFlow.",
      "action": "move-key",
      "parameters": {
        "source_path": "/problem/primary_equation/solver",
        "destination_path": "/problem/primary_equation/nonlinear_solver/linear_solver",
        "create_path":true  
      }
    },
    {
      "action": "move-key",
      "parameters": {
        "source_path": "/problem/primary_equation/solver",
        "destination_path": "/problem/primary_equation/nonlinear_solver/linear_solver",
        "create_path":true  
      }
    },
    {
      "NAME": "Set storativity for Unsteady_LMH.",
      "action": "add-key",
      "parameters": {
        "path": "/problem/primary_equation",        
        "key": "aux_region_key",
        "value": "ALL",
        "path-type-filter" : "Unsteady_LMH",
        "path-type-filter-path" : "/problem/primary_equation"
      }
    },      
    {
      "NAME": "Set storativity for Unsteady_LMH.",
      "action": "add-key",
      "parameters": {
        "path": "/problem/primary_equation",        
        "key": "aux_storativity_key",
        "value": "1.0",
        "path-type-filter" : "Unsteady_LMH",
        "path-type-filter-path" : "/problem/primary_equation"
      }
    },      
    {
      "NAME": "Set storativity for Unsteady_LMH.",
      "action" : "move-key",
      "parameters": {
        "source_path":"/problem/primary_equation/aux_region_key",
        "destination_path":"/problem/primary_equation/_input_fields/0/region",
        "path-type-filter" : "Unsteady_LMH",
        "path-type-filter-path" : "/problem/primary_equation",
        "create_path":true
      }  
    },      
    {
      "NAME": "Set storativity for Unsteady_LMH.",
      "action" : "move-key",
      "parameters": {
        "source_path":"/problem/primary_equation/aux_storativity_key",
        "destination_path":"/problem/primary_equation/_input_fields/0/storativity",
        "path-type-filter" : "Unsteady_LMH",
        "path-type-filter-path" : "/problem/primary_equation",
        "create_path":true
      }  
    },
    {
      "NAME": "Set storativity for Unsteady_LMH.",
      "action" : "merge-arrays",
      "parameters": {
        "source_path":"/problem/primary_equation/_input_fields",
        "addition_path":"/problem/primary_equation/input_fields",
        "destination_path":"/problem/primary_equation/input_fields",
        "path-type-filter" : "Unsteady_LMH",
        "path-type-filter-path" : "/problem/primary_equation"     
      }
    },
    {
      "NAME": "Set storativity for Unsteady_LMH.",
      "action" : "move-key-forward",
      "parameters": {
        "path":"/problem/primary_equation/input_fields",
        "path-type-filter" : "Unsteady_LMH",
        "path-type-filter-path" : "/problem/primary_equation"
      }
    },
    {
      "NAME": "Set storativity for Unsteady_MH.",
      "action": "add-key",
      "parameters": {
        "path": "/problem/primary_equation",        
        "key": "aux_region_key",
        "value": "ALL",
        "path-type-filter" : "Unsteady_MH",
        "path-type-filter-path" : "/problem/primary_equation"
      }
    },      
    {
      "NAME": "Set storativity for Unsteady_MH.",
      "action": "add-key",
      "parameters": {
        "path": "/problem/primary_equation",        
        "key": "aux_storativity_key",
        "value": "1.0",
        "path-type-filter" : "Unsteady_MH",
        "path-type-filter-path" : "/problem/primary_equation"
      }
    },      
    {
      "NAME": "Set storativity for Unsteady_MH.",
      "action" : "move-key",
      "parameters": {
        "source_path":"/problem/primary_equation/aux_region_key",
        "destination_path":"/problem/primary_equation/_input_fields/0/region",
        "path-type-filter" : "Unsteady_MH",
        "path-type-filter-path" : "/problem/primary_equation",
        "create_path":true  
            
      }
    },      
    {
      "NAME": "Set storativity for Unsteady_MH.",
      "action" : "move-key",
      "parameters": {
        "source_path":"/problem/primary_equation/aux_storativity_key",
        "destination_path":"/problem/primary_equation/_input_fields/0/storativity",
        "path-type-filter" : "Unsteady_MH",
        "path-type-filter-path" : "/problem/primary_equation",
        "create_path":true  
            
      }
    },      
    {
      "NAME": "Set storativity for Unsteady_MH.",
      "action" : "merge-arrays",
      "parameters": {
        "source_path":"/problem/primary_equation/_input_fields",
        "addition_path":"/problem/primary_equation/input_fields",
        "destination_path":"/problem/primary_equation/input_fields",
        "path-type-filter" : "Unsteady_MH",
        "path-type-filter-path" : "/problem/primary_equation"
      }
    },            
    {
      "NAME": "Set storativity for Unsteady_MH.",
      "action" : "move-key-forward",
      "parameters": {
        "path":"/problem/primary_equation/input_fields",
        "path-type-filter" : "Unsteady_MH",
        "path-type-filter-path" : "/problem/primary_equation"
      }
    },
    {
      "NAME": "Use time aware Darcy_MH instead of Steady.",
      "action": "rename-type",
      "parameters": {
        "path": "/problem/primary_equation",
        "old_name": "Steady_MH",
        "new_name": "Flow_Darcy_MH"
      }
    },
        {
      "NAME": "Use time aware Darcy_MH instead of Steady.",
      "action": "rename-type",
      "parameters": {
        "path": "/problem/primary_equation",
        "old_name": "Unsteady_MH",
        "new_name": "Flow_Darcy_MH"
      }
    },
    {
      "NAME": "Use time aware Darcy_MH instead of Steady.",
      "action": "rename-type",
      "parameters": {
        "path": "/problem/primary_equation",
        "old_name": "Unsteady_LMH",
        "new_name": "Flow_Richards_LMH"
      }
    },
    {
      "NAME": "Rename sequential coupling",
      "action": "rename-type",
      "parameters": {
        "path": "/problem",
        "old_name": "SequentialCoupling",
        "new_name": "Coupling_Sequential"
      }
    },
    {
      "NAME": "Rename sequential coupling keys",
      "action": "move-key",
      "parameters": {
        "source_path": "/problem/primary_equation",
        "destination_path": "/problem/flow_equation"
      }
    },
    {
      "NAME": "Rename sequential coupling keys",
      "action": "move-key",
      "parameters": {
        "source_path": "/problem/secondary_equation",
        "destination_path": "/problem/solute_equation",
        "type-filter": "Coupling_OperatorSplitting"
      }
    },
    {
      "NAME": "Rename sequential coupling keys",
      "action": "move-key",
      "parameters": {
        "source_path": "/problem/secondary_equation",
        "destination_path": "/problem/heat_equation",
        "type-filter": "Heat_AdvectionDiffusion_DG"
      }
    },

    {
      "NAME": "Make output-specific key.",
      "action": "add-key",
      "parameters": {
        "path": "/problem/flow_equation",
        "key": "output_specific"
      }
    },
    {
      "NAME": "Move to output specific.",
      "action": "move-key",
      "parameters": {
        "source_path": "/problem/flow_equation/output/raw_flow_output",
        "destination_path": "/problem/flow_equation/output_specific/raw_flow_output",
        "create_path":true
      }
    },
    {
      "NAME": "Move to output specific.",
      "action": "move-key",
      "parameters": {
        "source_path": "/problem/flow_equation/output/compute_errors",
        "destination_path": "/problem/flow_equation/output_specific/compute_errors",
        "create_path":true
      }
    },
    {
      "NAME": "Move DarcyFlow output_stream.",
      "action": "move-key",
      "parameters": {
        "source_path": "/problem/flow_equation/output/output_stream",
        "destination_path": "/problem/flow_equation/output_stream",
        "create_path":true
      }
    },
    {
      "NAME": "Rename DarcyFlow output_fields.",
      "action": "move-key",
      "parameters": {
        "source_path": "/problem/flow_equation/output/output_fields",
        "destination_path": "/problem/flow_equation/output/fields"
      }
    },
    {
      "NAME": "Move time step for DarcyFlow output stream.",
      "action": "move-key",
      "parameters": {
        "source_path": "/problem/flow_equation/output_stream/time_step",
        "destination_path": "/problem/flow_equation/output_stream/times/0/step",
        "create_path":true
      }
    },
    {
      "NAME": "Move time_list for Darcy.",
      "action" : "merge-arrays",
      "parameters": {
        "source_path": "/problem/flow_equation/output_stream/time_list",
        "addition_path": "/problem/flow_equation/output_stream/times",        
        "destination_path": "/problem/flow_equation/output_stream/times"
      }
    },    
    {
      "NAME": "Move time step for DarcyFlow output stream.",
      "action": "move-key",
      "parameters": {
        "source_path": "/problem/flow_equation/output_stream/add_input_times",
        "destination_path": "/problem/flow_equation/output/add_input_times",
        "create_path":true
      }
    },




    {
      "NAME": "Make output_fields in transport.",
      "action": "move-key",
      "parameters": {
        "source_path": "/problem/solute_equation/transport/output_fields",
        "destination_path": "/problem/solute_equation/transport/output/fields",
        "create_path":true
      }
    },
    {
      "NAME": "Make output_fields in dual porosity.",
      "action": "move-key",
      "parameters": {
        "source_path": "/problem/solute_equation/reaction_term/output_fields",
        "destination_path": "/problem/solute_equation/reaction_term/output/fields",
        "create_path":true
      }
    },
    {
      "NAME": "Make output_fields in mobile reaction.",
      "action": "move-key",
      "parameters": {
        "source_path": "/problem/solute_equation/reaction_term/reaction_mobile/output_fields",
        "destination_path": "/problem/solute_equation/reaction_term/reaction_mobile/output/fields",
        "create_path":true
      }
    },
    {
      "NAME": "Make output_fields in immobile reaction.",
      "action": "move-key",
      "parameters": {
        "source_path": "/problem/solute_equation/reaction_term/reaction_immobile/output_fields",
        "destination_path": "/problem/solute_equation/reaction_term/reaction_immobile/output/fields",
        "create_path":true
      }
    },
    {
      "NAME": "Make time step for transport output stream.",
      "action": "move-key",
      "parameters": {
        "source_path": "/problem/solute_equation/output_stream/time_step",
        "destination_path": "/problem/solute_equation/output_stream/times/0/step",
        "create_path":true
      }
    },
    {
      "NAME": "Move time_list for transport.",
      "action" : "merge-arrays",
      "parameters": {
        "source_path": "/problem/solute_equation/output_stream/time_list",
        "addition_path": "/problem/solute_equation/output_stream/times",        
        "destination_path": "/problem/solute_equation/output_stream/times"
      }
    },    
    {
      "NAME": "Move add_input_times for transport.",
      "action": "move-key",
      "parameters": {
        "source_path": "/problem/solute_equation/output_stream/add_input_times",
        "destination_path": "/problem/solute_equation/output/add_input_times",
        "keep_source":true,
        "create_path":true
      }
    },      
    {
      "NAME": "Move add_input_times for transport.",
      "action": "move-key",
      "parameters": {
        "source_path": "/problem/solute_equation/output_stream/add_input_times",
        "destination_path": "/problem/solute_equation/reaction_term/output/add_input_times",
        "keep_source":true,
        "create_path":true
      }
    },
    {
      "NAME": "Move add_input_times for transport.",
      "action": "move-key",
      "parameters": {
        "source_path": "/problem/solute_equation/output_stream/add_input_times",
        "destination_path": "/problem/solute_equation/reaction_term/reaction_mobile/output/add_input_times",
        "keep_source":true,
        "create_path":true
      }
    },
    {
      "NAME": "Move add_input_times for transport.",
      "action": "move-key",
      "parameters": {
        "source_path": "/problem/solute_equation/output_stream/add_input_times",
        "destination_path": "/problem/solute_equation/reaction_term/reaction_immobile/output/add_input_times",        
        "create_path":true
      }
    },    
      




    {
      "NAME": "Rename Heat output_fields.",
      "action": "move-key",
      "parameters": {
        "source_path": "/problem/heat_equation/output/output_fields",
        "destination_path": "/problem/heat_equation/output/fields",
        "create_path":true
      }
    },
    {
      "NAME": "Move time step for Heat output stream.",
      "action": "move-key",
      "parameters": {
        "source_path": "/problem/heat_equation/output_stream/time_step",
        "destination_path": "/problem/heat_equation/output_stream/times/0/step",
        "create_path":true
      }
    },
    {
      "NAME": "Move time_list for Heat.",
      "action" : "merge-arrays",
      "parameters": {
        "source_path": "/problem/heat_equation/output_stream/time_list",
        "addition_path": "/problem/heat_equation/output_stream/times",        
        "destination_path": "/problem/heat_equation/output_stream/times"
      }
    },    
    {
      "NAME": "Move time step for Heat output stream.",
      "action": "move-key",
      "parameters": {
        "source_path": "/problem/heat_equation/output_stream/add_input_times",
        "destination_path": "/problem/heat_equation/output/add_input_times",
        "create_path":true
      }
    },    
    {
        "NAME": "Change balance:true",
        "action": "change-value",
        "parameters":{
          "path" : "/problem/flow_equation/balance",
          "old_value" : "true",
          "new_value" : "{}"
        }
    },
    {
        "NAME": "Change balance:true",
        "action": "change-value",
        "parameters":{
          "path" : "/problem/solute_equation/balance",
          "old_value" : "true",
          "new_value" : "{}"
        }
    },    
    {
        "NAME": "Change balance:true",
        "action": "change-value",
        "parameters":{
          "path" : "/problem/heat_equation/balance",
          "old_value" : "true",
          "new_value" : "{}"
        }
    },    
    {
        "NAME": "Change balance:true",
        "action": "change-value",
        "parameters":{
          "path" : "/problem/flow_equation/balance",
          "old_value" : "false",
          "new_value" : "{add_output_times: false}"
        }
    },
    {
        "NAME": "Change balance:true",
        "action": "change-value",
        "parameters":{
          "path" : "/problem/solute_equation/balance",
          "old_value" : "false",
          "new_value" : "{add_output_times: false}"
        }
    },    
    {
        "NAME": "Change balance:true",
        "action": "change-value",
        "parameters":{
          "path" : "/problem/heat_equation/balance",
          "old_value" : "false",
          "new_value" : "{add_output_times: false}"
        }
    },
    {
        "NAME": "Change BOUNDARY to .BOUNDARY",
        "action": "change-value",
        "parameters":{
          "path" : "/problem/*/input_fields/*/region",
          "old_value" : "BOUNDARY",
          "new_value" : ".BOUNDARY"
        }
    },
    {
        "NAME": "Change BOUNDARY to .BOUNDARY",
        "action": "change-value",
        "parameters":{
          "path" : "/problem/solute_equation/*/input_fields/*/region",
          "old_value" : "BOUNDARY",
          "new_value" : ".BOUNDARY"
        }
    },
    {
        "NAME": "Change BOUNDARY to .BOUNDARY",
        "action": "change-value",
        "parameters":{
          "path" : "/problem/solute_equation/reaction_term/*/input_fields/*/region",
          "old_value" : "BOUNDARY",
          "new_value" : ".BOUNDARY"
        }
    },
    {
        "NAME": "Change BOUNDARY to .BOUNDARY, hack to deal with substitution matching the substrings",
        "action": "change-value",
        "parameters":{
          "path" : "/problem/*/input_fields/*/region",
          "old_value" : "IMPLICIT .BOUNDARY",
          "new_value" : ".IMPLICIT_BOUNDARY"
        }
    },
    {
        "NAME": "Change BOUNDARY to .BOUNDARY",
        "action": "change-value",
        "parameters":{
          "path" : "/problem/solute_equation/*/input_fields/*/region",
          "old_value" : "IMPLICIT .BOUNDARY",
          "new_value" : ".IMPLICIT_BOUNDARY"
        }
    },
    {
        "NAME": "Change BOUNDARY to .BOUNDARY",
        "action": "change-value",
        "parameters":{
          "path" : "/problem/solute_equation/reaction_term/*/input_fields/*/region",
          "old_value" : "IMPLICIT .BOUNDARY",
          "new_value" : ".IMPLICIT_BOUNDARY"
        }
    }
            
             
  ]
}

'''

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--from-version", help="Version of the input file.")
    parser.add_argument("-t", "--to-version", help="Version of the output.")
    parser.add_argument("-r", "--reverse", action='store_true', help="Perform reversed conversion. Input file is in 'to-version'.")
    parser.add_argument('in_file', help="Input YAML (or CON) file(s). Wildcards accepted.")

    args = parser.parse_args()


    path_wild = args.in_file.split('/')
    if path_wild[0]:
        # relative path
        dirs = ["."]
    else:
        # absloute
        dirs = ["/"]

    files=[]
    for wild_name in path_wild:
        # separate dirs and files
        paths=[]
        for dir in dirs:
            paths += [ dir + "/" + f for f in os.listdir(dir) if fnmatch.fnmatch(f, wild_name) ]
        dirs=[]
        for path in paths:
            if os.path.isdir(path):
                dirs.append(path)
            elif os.path.isfile(path):
                files.append(path)
            else:
                assert False, "Path neither dir nor file."

    print(files)
    trees = []
    for fname in files:
        with open(fname, "r") as f:
            trees.append( ( fname, yml.load(f)) )

    changes.apply_changes(trees, args.from_version, args.to_version, reversed=args.reverse)

    for fname, tree in trees:
        base=os.path.splitext(fname)[0]
        out_fname = base + ".new.yaml"
        with open(out_fname, "w") as f:
            yml.dump(tree, f)