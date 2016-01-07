"""
File for transformation of yaml files to another version

Wildchars:
    In defined actions (if is not writen in next text else) is possible use simple
    wildchars. For key with name path and source_path:
    - * for one level of path
    - ** for one or more level of path
    For key with name destination_path, addition_path and set_type_path is 
    possible use replacement $1,$2,... .  This signify part of path in wildchars 
    on position 1,2,....
    Wildchars have three restrictions:
        - ** can't be last part of path (too many result cases)
        - **/* sequence is forbiden (too many result cases)
        - **/** sequence is forbiden (can't determine replacements)
        
Filter:
    - type-filter is a optional parameter for a key in possition path or source_path.
      If is the parameter set, a opperation is processed only for the key of set type.
    - parent-type-filter is a optional filter parameter similar as type-filter but for
      parent of set key      
    
Actions:
    - move-key - Move value of set key from source_path to destination_path.
      Not existing directory of destination_path cause exception. If both path
      is same, only key is different, key is renamed. Destination_path placed in
      source_path cause exception. If source_path not exist action si skipped.
      Refference or anchors is moved and if its relative possition is changed,
      result must be fixed by user. In this action can be used parameters
      set_type_path and new_type for replacing tag new_type in path "set_type_path"
    - delete-key - Delete key on set path. If key contains anchor or refference,
      transporter try resolve reference. If path not exist, or contains some children
      action si skipped. If parameter deep is set to true, key is delete with containing 
      children.
    - rename-type - Change tag on set path from old_name to new_name.
    - change-value - Change value on set path from old_value to new_value.
      If node in specified path is not scallar type node, exception is raise.
    - scale-value - multiply value in set path by scale. If node in specified 
      path is not scallar type node with numeric value, exception is raise.
    - merge-arrays - Add array in addition_path after array in source path. If
      destination_path is defined, node is moved from source to this path.
    - move-key-forward - Move node in specified path before all siblings nodes.
      (this action is use by import script for moving node with anchor after
      refference)
    - add-key - add set key to path. Key parent should be existing abstract record.
      if optional parameters value or type is set, corresponding key parameter is 
      added.

Description:
    Transformator check se json transformation file. If this file is in bad format,
    cause exception. After checking is processet action in set order. Repetitional
    processing this program can cause exception or corrupt data.
    
Example::

    {
        "name": "Example file",
        "description": "Template for creating transformation file",
        "old_format": "",
        "new_format": "",
        "actions": [
            {
                "action": "move-key",
                "parameters":{
                    "source_path":"/problem/primary_equation/input_fields",
                    "destination_path":"/problem/mesh/balance"
                }
            },
            {
                "action": "delete-key",
                "parameters": {
                    "path": "/problem/mesh/mesh_file"
                }
            },
            {
                "action": "rename-type",
                "parameters": {
                    "path": "/problem/primary_equation",
                    "old_name": "Steady_MH",
                    "new_name": "SteadyMH"
                }
            }
        ]
    }

Wildchars Example::

    { 
      "action": "move-key",
      "parameters": {
        "source_path":"/**/input_fields/*/r_set",
        "destination_path":"/$1/input_fields/$2/region"
      }
    }

"""

# pylint: disable=invalid-name

import json
import re
import copy

from .loader import Loader
from ..data_node import DataNode
from helpers import NotificationHandler
from helpers.subyaml.structure_changer import StructureChanger
from data.autoconversion import autoconvert
from data.format import get_root_input_type_from_json

class Transformator:
    """Transform yaml file to new version"""
    __actions__=["delete-key", "move-key", "rename-type", "move-key-forward",
                          "change-value", "merge-arrays", "add-key", "scale-value"]
    __source_paths__ = {"delete-key":"path","move-key-forward":"path", "move-key":"source_path", 
                                      "rename-type":"path", "change-value":"path", "merge-arrays":"source_path", 
                                      "add-key":"path", "scale-value":"path"}
    __destination_paths__ = {"delete-key":None,"move-key-forward":None, "move-key":["destination_path"],
                                          "rename-type":None, "change-value":None, "merge-arrays":["destination_path",
                                          "addition_path"], "add-key":None, "scale-value":None}
    
    def __init__(self, transform_file, data=None):
        """init"""
        if transform_file is not None:
            self._transformation = json.loads(transform_file)
        else:
            self._transformation = data
        "parsed json transformation file"
        # check
        if not 'actions' in  self._transformation:
            raise TransformationFileFormatError("Array of actions is required")
        i = 1
        for action in self._transformation['actions']:
            self._check_parameter("action", action, action['action'], i)
            if not action['action'] in Transformator.__actions__:
                raise TransformationFileFormatError(
                    " Action '" + action['action'] + "' is not specified")
            self._check_parameter("parameters", action, action['action'], i)
            if action['action'] == "delete-key":
                self._check_parameter("path", action['parameters'], action['action'], i)
            elif action['action'] == "move-key-forward":
                self._check_parameter("path", action['parameters'], action['action'], i)
            elif action['action'] == "move-key":
                self._check_parameter("destination_path", action['parameters'], action['action'], i)
                self._check_parameter("source_path", action['parameters'], action['action'], i)
                if 'create_path' not in action['parameters']:
                    action['parameters']['create_path'] = False
                if "set_type_path" in action['parameters'] or "new_type" in action['parameters']:
                    self._check_parameter("set_type_path", action['parameters'], action['action'], i)
                    self._check_parameter("new_type", action['parameters'], action['action'], i)                
            elif action['action'] == "rename-type":
                self._check_parameter("path", action['parameters'], action['action'], i)
                self._check_parameter("new_name", action['parameters'], action['action'], i)
                self._check_parameter("old_name", action['parameters'], action['action'], i)
            elif action['action'] == "change-value":
                self._check_parameter("path", action['parameters'], action['action'], i)
                self._check_parameter("new_value", action['parameters'], action['action'], i)
                self._check_parameter("old_value", action['parameters'], action['action'], i)
            elif action['action'] == "merge-arrays":
                self._check_parameter("source_path", action['parameters'], action['action'], i)
                self._check_parameter("addition_path", action['parameters'], action['action'], i)
            elif action['action'] == "add-key":
                self._check_parameter("path", action['parameters'], action['action'], i) 
                self._check_parameter("key", action['parameters'], action['action'], i)
            elif action['action'] == "scale-value":
                self._check_parameter("path", action['parameters'], action['action'], i)
                self._check_parameter("scale", action['parameters'], action['action'], i)
                try:
                    float(action['parameters']['scale'])
                except ValueError:
                    raise TransformationFileFormatError(
                        "Type of parameter scale for action type scale-value is not " +
                        "numeric (value: " + action['parameters']['scale'] + ")")
            i += 1

    def _check_parameter(self, name, dict_, act_type, line):
        if name not in dict_:
            raise TransformationFileFormatError(
                "Parameter " + name + " for action type '" + act_type + "' si required" +
                " (action: " + str(line) + ")")
        if dict_[name] is None or dict_[name] == "":
            raise TransformationFileFormatError(
                "Parameter " + name + " for action type '" + act_type + "' cannot be empty" +
                " (action: " + str(line) + ")")

    @property
    def old_version(self):
        """return version of yaml document befor transformation"""
        if 'old_format' in self._transformation:
            return self._transformation['old_format']
        return ""

    @property
    def new_version(self):
        """return version of yaml document after transformation"""
        if 'new_format' in self._transformation:
            return self._transformation['new_format']
        return ""

    @property
    def description(self):
        """return transformation description"""
        if 'description' in self._transformation:
            return self._transformation['description']
        return ""

    @property
    def name(self):
        """return transformation name"""
        if 'name' in self._transformation:
            return self._transformation['name']
        return ""

    def transform(self, yaml, cfg):
        """transform yaml file"""
        # TODO: cls.root = autoconvert(cls.root, cls.root_input_type)
        notification_handler = NotificationHandler()
        loader = Loader(notification_handler)
        root = loader.load(yaml)
        lines = yaml.splitlines(False)
        changes = False
        text = cfg.get_curr_format_text()
        root_input_type = get_root_input_type_from_json(text)
        for aaction in self._transformation['actions']:
            if changes:
                root, lines = self.refresh(root_input_type, yaml, notification_handler, loader)
                changes=False
            for action in self._replace_wildchars(root, aaction):
                if changes:
                    root, lines = self.refresh(root_input_type, yaml, notification_handler, loader)
                if action['action'] == "delete-key":
                    changes = self._delete_key(root, lines, action)
                elif action['action'] == "move-key":
                    changes = self._move_key(root, lines, action)
                    if  changes and "set_type_path" in action['parameters']:
                        yaml = "\n".join(lines)
                        root, lines = self.refresh(root_input_type, yaml, notification_handler, loader)
                        changes = False                        
                        try:
                            node = root.get_node_at_path(action['parameters']['set_type_path']) 
                            if node.type is not None:
                                StructureChanger.change_tag(lines, node, node.type, action['parameters']['new_type'])
                            else:
                                StructureChanger.add_tag(lines, node, action['parameters']['new_type'])                                
                            changes = True                                
                        except:
                            pass
                elif action['action'] == "rename-type":
                    changes = self._rename_type(root, lines, action)
                elif action['action'] == "move-key-forward":
                    changes = self._move_key_forward(root, lines, action)
                elif action['action'] == "change-value":
                    changes = self._change_value(root, lines, action)
                elif action['action'] == "merge-arrays":
                    changes = self._add_array(root, lines, action)
                    if 'destination_path' in action['parameters']:
                        if changes: 
                            yaml = "\n".join(lines)
                            root, lines = self.refresh(root_input_type, yaml, notification_handler, loader)
                        changes = self._move_key(root, lines, action)
                elif action['action'] == "add-key":
                    changes = self._add_key(root, lines, action)
                elif action['action'] == "scale-value":
                    changes = self._scale_value(root, lines, action)
                if changes:
                    yaml = "\n".join(lines)
        yaml = "\n".join(lines)
        return yaml
        
    def refresh(self, root_input_type, yaml, notification_handler, loader):
        notification_handler.clear()
        root = loader.load(yaml)
        autoconvert(root, root_input_type)
        lines = yaml.splitlines(False)        
        return root,  lines

    class _Replacement:
        """Replacement structure for path wildchars operation"""
        def __init__(self, replacement=None):
            self.path = ""
            """find path"""
            self.replacements = []
            """array of replacements"""
            self.orig_path = None
            """origination path for error message"""
            if replacement is not None:
                self.path = copy.deepcopy(replacement.path)
                self.replacements = copy.deepcopy(replacement.replacements)
                self.orig_path = copy.deepcopy(replacement.orig_path)

    def _replace_wildchars(self, root, action):
        """If path have some wildchards replace it real path"""        
        if Transformator.__source_paths__[action['action']] is None:
            return[action]
        path_parameter = Transformator.__source_paths__[action['action']]
        path = action['parameters'][path_parameter]
        if '*' not in path:
            if not self._filtered(path, action, root):
                return []
            return [action]
        res = []
        spath = path.split('/')
        if spath[0] == "":
            spath = spath[1:]
        if spath[len(spath)-1]=='**':
            raise TransformationFileFormatError(
                "Wildcard '**' can't be in end of path(" + action['parameters'][path_parameter] + ")")
        replacement = Transformator._Replacement()
        replacement.orig_path = path
        replacements = self._search_node(spath, root, replacement, action['parameters'][path_parameter])
        for replacement in replacements:
            if not self._filtered(replacement.path, action, root):
                continue
            new_action = copy.deepcopy(action)
            new_action['parameters'][path_parameter] = replacement.path 
            new_action['parameters']['orig_' + path_parameter] = replacement.orig_path
            if Transformator.__destination_paths__[action['action']] is not None:
                for parameter in Transformator.__destination_paths__[action['action']]:
                    if parameter in action['parameters']:
                        for i in range(0, len(replacement.replacements)):
                            if i == 0:
                                new_action['parameters']['orig_' + parameter] = new_action['parameters'][parameter]
                            new_action['parameters'][parameter] = \
                                new_action['parameters'][parameter].replace( \
                                    '$'+str(i+1), replacement.replacements[i]) 
                            if new_action['parameters'][parameter][0] != "/":
                                new_action['parameters'][parameter] = "/"+new_action['parameters'][parameter]
            res.append(new_action)
        if action['action'] == "delete-key": 
            # keys in array must be delete from end
            res.reverse()
        return res

    def _filtered(self,  path, action, root, parent=False):
        """Return if type of set record have requaired type"""
        if 'type-filter' in action['parameters']:
            try:
                node = root.get_node_at_path(path)
                if node.type.value != action['parameters']['type-filter']:
                    return False
            except:
                return False
        if 'parent-type-filter' in action['parameters']:
            try:
                node = root.get_node_at_path(path).parent
                if node.type.value != action['parameters']['parent-type-filter']:
                    return False
            except:
                return False
        return True
        

    def _search_node(self, spath, node, replacement,  orig_path, deep=False):
        """
        search all paths complying set spath pattern in set node, 
        and return theirs array. This function is use recusively.
        """
        if node is None:
            return []
        if len(spath) == 0:
            return [replacement]
        res = []
        if spath[0]=='*':
            for child in node.children_keys:
                new_replacement = Transformator._Replacement(replacement)
                new_replacement.path += "/" + child
                new_replacement.replacements.append(child)
                if len(spath) == 1:
                    res.append(new_replacement)
                else:
                    replacements = self._search_node(spath[1:],  node.get_child(child), new_replacement, orig_path)                
                    res.extend(replacements)
        elif spath[0]=='**':
            if len(spath)>1:
                if spath[1]=='**':
                    raise TransformationFileFormatError(
                        "Sequence '**/**' is forbiden(" + orig_path + ")")
                if spath[1]=='*':
                    raise TransformationFileFormatError(
                        "Sequence '**/*' is forbiden(" + orig_path + ")")
            for child in node.children_keys:
                new_replacement = Transformator._Replacement(replacement)
                new_replacement.path += "/" + child
                if deep:
                    new_replacement.replacements[len(new_replacement.replacements)-1] += \
                        "/" + child
                else:
                    new_replacement.replacements.append(child)
                if len(spath) == 1:
                    res.append(new_replacement)
                # first try find paths in next structure
                replacements = self._search_node(spath, node.get_child(child), new_replacement, orig_path, True)                
                res.extend(replacements)
                if deep:
                    # second try cancel ** find
                    if spath[1] == child:
                        new_replacement2 = Transformator._Replacement(replacement)
                        new_replacement2.path += "/" + child
                        if len(spath)==2:
                            res.append(new_replacement2)
                        else:
                            replacements = self._search_node(spath[2:],  node.get_child(child), new_replacement2, orig_path)                
                            res.extend(replacements)
        else:
            for child in node.children_keys:
                if child == spath[0]:
                    new_replacement = Transformator._Replacement(replacement)
                    new_replacement.path += "/" + child
                    if len(spath) == 1:
                        res.append(new_replacement)
                    else:
                        replacements = self._search_node(spath[1:], node.get_child(child), new_replacement, orig_path)                
                        res.extend(replacements)
        return res

    def _find_all_ref(self, node, refs):
        """find all references"""
        for child in node.children:
            if child.anchor is not None and child.ref is not None:
                if child.anchor.value not in refs:
                    refs[child.anchor.value] = {}
                    refs[child.anchor.value]['anchor'] = child.ref
                    refs[child.anchor.value]['ref'] = [child]
                else:
                    refs[child.anchor.value]['ref'].append(child)
            if not child.implementation == DataNode.Implementation.scalar:
                self._find_all_ref(child, refs)

    def _move_key_forward(self, root, lines, action):
        """Move key forward"""
        try:
            parent = re.search(r'^(.*)/([^/]*)$', action['parameters']['path'])
            node = root.get_node_at_path(action['parameters']['path'])
            if parent is None:
                raise TransformationFileFormatError(
                    "Cannot find parent path for path (" + self._get_paths_str(action, 'path') + ")")
            is_root = False
            if len(parent.group(1)) == 0:
                parent_node = root
                is_root = True
            else:
                parent_node = root.get_node_at_path(parent.group(1))
        except:
            return False
        l1, c1, l2, c2 = StructureChanger.node_pos(node)
        pl1, pc1, pl2, pc2 = parent_node.span.start.line-1, parent_node.span.start.column-1, \
            parent_node.span.end.line-1, parent_node.span.end.column-1
        if parent_node.implementation != DataNode.Implementation.mapping:
            raise TransformationFileFormatError(
                "Parent of path (" + self._get_paths_str(action, 'path') + ") must be abstract record")
        if is_root:
            intendation1 = 0
            pl1 = 0
            pc1 = 0
        else:
            intendation1 = re.search(r'^(\s*)(\S.*)$', lines[pl1])
            intendation1 = len(intendation1.group(1)) + 2
        l1, c1, l2, c2 = StructureChanger._add_comments(lines, l1, c1, l2, c2)
        add = StructureChanger.copy_structure(lines, l1, c1, l2, c2, intendation1)
        pl1, pc1 = StructureChanger.skip_tar(lines, pl1, pc1, pl2, pc2)
        action['parameters']['deep'] = True
        self._delete_key(root, lines, action)
        StructureChanger.paste_structure(lines, pl1, add, pc1 != 0)
        return True
        
    def _get_paths_str(self, action, key):
        """return error string with both paths"""
        str = action['parameters'][key]
        if ('orig_' + key) in action['parameters'] and \
            action['parameters']['orig_' + key] is not None:
            str += " ~" + action['parameters']['orig_' + key] + "~"
        return str

    def _delete_key(self, root, lines, action):
        """Delete key transformation"""
        try:
            node = root.get_node_at_path(action['parameters']['path'])
        except:
            return False
        if "deep" not in action['parameters'] or not  action['parameters']['deep']:
            if len(node.children_keys)>0:
                return False            
        l1, c1, l2, c2 = StructureChanger.node_pos(node)
        anchors = []
        for i in range(l1, l2+1):
            if l1 == l2:
                text = lines[l1][c1:c2]
            elif i == l1:
                text = lines[l1][c1:]
            elif i == l2:
                text = lines[l2][:c2]
            else:
                text = lines[i]
            anchor = re.search(r'\s+&\s*(\S+)\s+', text)
            if anchor is None:
                anchor = re.search(r'\s+&\s*(\S+)$', text)
            if anchor is not None:
                anchors.append(anchor.group(1))
        if len(anchors) > 0:
            self._rewrite_first_ref(root, lines, anchors)
        # try delete with comma after or prev
        nl2, nc2 = self._next_spread(lines, l2, c2)
        if nl2 == l2 and nc2 == c2:
            l1, c1 = self._prev_spred(lines, l1, c1)
        else:
            l2 = nl2
            c2 = nc2
        # try exclude comments
        l1, c1, l2, c2 = StructureChanger.leave_comments(lines, l1, c1, l2, c2) 
        l1, c1, l2, c2 = StructureChanger.add_delete_item_chars(lines, l1, c1, l2, c2)
        StructureChanger.delete_structure(lines, l1, c1, l2, c2)
        return True

    def _delete_value(self, root, lines, node):
        """Delete value and return start possition deleted value"""
        l1, c1, l2, c2 = node.span.start.line-1, node.span.start.column-1, \
            node.span.end.line-1, node.span.end.column-1
        l1, c1 = StructureChanger.skip_tar(lines, l1, c1, l2, c2)
        StructureChanger.delete_structure(lines, l1, c1, l2, c2)
        return l1, c1

    def _rewrite_first_ref(self, root, lines, anchors):
        """copy value from anchor node to first refference"""
        refs = {}
        self._find_all_ref(root, refs)
        refs_anchor = []
        refs_ref = []
        refs_remove = []
        for r in refs:
            if r in anchors and refs[r]['anchor'] is not None and len(refs[r]['ref']) > 0:
                ref_node = refs[r]['ref'][0]
                for ref in refs[r]['ref']:
                    if (ref_node.span.start.line > ref.span.start.line or
                            (ref_node.span.start.line == ref.span.start.line and
                             ref_node.span.start.column > ref.span.start.column)):
                        ref_node = ref
                added = False
                for i in range(0, len(refs_anchor)):
                    if (ref_node.span.start.line > refs_ref[i].span.start.line or
                            (ref_node.span.start.line == refs_ref[i].span.start.line and
                             ref_node.span.start.column > refs_ref[i].span.start.column)):
                        refs_anchor.insert(i, refs[r]['anchor'])
                        refs_ref.insert(i, ref_node)
                        refs_remove.insert(i, len(refs[r]['ref']) == 1)
                        added = True
                        break
                if not added:
                    refs_anchor.append(refs[r]['anchor'])
                    refs_ref.append(ref_node)
                    refs_remove.append(len(refs[r]['ref']) == 1)
        for i in range(0, len(refs_anchor)):
            self._copy_tree_from_anchor(lines, refs_anchor[i], refs_ref[i], refs_remove[i])

    def _copy_tree_from_anchor(self, lines, anchor_node, ref_node, remove_anchor=False):
        """Copy tree structure from anchor node to node with reference"""
        l1, c1, l2, c2 = anchor_node.span.start.line-1, anchor_node.span.start.column-1, \
            anchor_node.span.end.line-1, anchor_node.span.end.column-1
        l1, c1 = StructureChanger.skip_tar(lines, l1, c1, l2, c2)
        hlpl1, hlpc1, l2, c2 = StructureChanger._add_comments(lines, l1, c1, l2, c2)
        dl1, dc1, dl2, dc2 = StructureChanger.node_pos(ref_node)
        intend = re.search(r'^(\s*)(\S.*)$', lines[dl1])
        intend = len(intend.group(1)) + 2
        add = StructureChanger.copy_structure(lines, l1, c1, l2, c2, intend)
        while dl1 <= dl2:
            ref = re.search(r'^(.*\*' + anchor_node.anchor.value + r')(.*)$', lines[dl1])
            if ref is not None:
                anchor = "&" + anchor_node.anchor.value
                if remove_anchor:
                    anchor = ""
                lines[dl1] = re.sub(r'\*' + anchor_node.anchor.value + r"\s+", anchor + " ",
                                    lines[dl1])
                lines[dl1] = re.sub(r'\*' + anchor_node.anchor.value + r"$", anchor, lines[dl1])
                if not ref.group(2).isspace() and len(ref.group(2)) > 0:
                    lines.insert(dl1+1, ref.group(2))
                    lines[dl1] = ref.group(1)
                StructureChanger.paste_structure(lines, dl1, add, True, True)                
                break
            dl1 += 1

    @staticmethod
    def _next_spread(lines, line, column):
        """if next no-empty char is comma, return new possition"""
        if line < len(lines):
            return line, column
        old_line = line
        old_column = column
        while True:
            column += 1
            while len(lines[line]) <= column:
                column = 0
                line += 1
                if line < len(lines):
                    break
            if not lines[line][column].isspace():
                if lines[line][column] == ",":
                    return line, column
                break
        return old_line, old_column

    @staticmethod
    def _prev_spred(lines, line, column):
        """if previous no-empty char is comma, return new position"""
        old_line = line
        old_column = column
        while True:
            column -= 1
            while 0 > column:
                line -= 1
                if line < 0:
                    break
                column = len(lines[line])-1
            if not lines[line][column].isspace():
                if lines[line][column] == ",":
                    return line, column
                break
        return old_line, old_column
  
    def _move_key(self, root, lines, action):
        """Move key transformation"""
        try:
            parent1 = root.get_node_at_path(action['parameters']['destination_path'])
            raise TransformationFileFormatError(
                "Destination path (" + self._get_paths_str(action, 'destination_path') + ") already exist")
        except:
            pass
        try:
            parent1 = re.search(r'^(.*)/([^/]*)$', action['parameters']['source_path'])
            node1 = root.get_node_at_path(action['parameters']['source_path'])
        except:
            return False
        try:
            parent2 = re.search(r'^(.*)/([^/]*)$', action['parameters']['destination_path'])
            new_node = parent2.group(2)
            node2 = root.get_node_at_path(parent2.group(1))
            node_struct = []
        except:
            if action['parameters']['create_path'] and parent2 is not None:
                node_struct,new_path = StructureChanger.copy_absent_path(root, 
                    lines, node1.parent, parent2.group(1))
                if len(node_struct) == 0:
                    raise TransformationFileFormatError(
                        "Can't constract destination path (" + 
                        self._get_paths_str(action, 'destination_path') + ")")
                try:                    
                    new_node = parent2.group(2)
                    parent2 = re.search(r'^(.*)/([^/]*)$', new_path + "/" + new_node)
                    node2 = root.get_node_at_path(parent2.group(1))                    
                except:    
                    raise TransformationFileFormatError(
                        "Constracted path error(" + self._get_paths_str(action, 'destination_path') +
                        ")")
            else:
                raise TransformationFileFormatError(
                    "Parent of destination path (" + self._get_paths_str(action, 'destination_path') +
                    ") must exist")
        sl1, sc1, sl2, sc2 = StructureChanger.node_pos(node1)
        dl1, dc1, dl2, dc2 = StructureChanger.node_pos(node2)
        if parent1.group(1) == parent2.group(1):
            # rename
            i = node1.key.span.start.line-1
            lines[i] = re.sub(parent1.group(2) + r"\s*:", parent2.group(2) + ":", lines[i])
            return True
        if not action['parameters']['create_path']:
            # check only existing path
            if node2.implementation != DataNode.Implementation.mapping:
                raise TransformationFileFormatError(
                    "Parent of destination path (" + self._get_paths_str(action, 'destination_path') +
                    ") must be abstract record")
        if node1.parent is None:
            intendation1 = 0
        else:
            intendation1 = re.search(r'^(\s*)(\S.*)$', lines[dl1])
            intendation1 = len(intendation1.group(1)) + 2
        
        sl1, sc1, sl2, sc2 = StructureChanger._add_comments(lines, sl1, sc1, sl2, sc2)
        add = StructureChanger.copy_structure(lines, sl1, sc1, sl2, sc2, intendation1)
        # rename key
        i = node1.key.span.start.line - sl1 - 1
        add[i] = re.sub(parent1.group(2) + r"\s*:", new_node + ":", add[i])        
        if action['parameters']['create_path'] and len(node_struct) > 0:
            add = StructureChanger.paste_absent_path(add, node_struct)        
        if sl1<dl1 and sl2>dl1:
            raise TransformationFileFormatError(
                "Destination block (" + self._get_paths_str(action, 'destination_path') +
                ") and source block (" + self._get_paths_str(action, 'source_path') +
                " is overlapped")
        if sl1 < dl2:
            # source before dest, first copy
            intendation2 = re.search(r'^(\s*)(\S.*)$', lines[dl2])
            StructureChanger.paste_structure(lines, dl2, add, len(intendation2.group(1)) < dc2)
            action['parameters']['path'] = action['parameters']['source_path']
            action['parameters']['deep'] = True
            self._delete_key(root, lines, action)
        else:
            # source after dest, first delete
            action['parameters']['path'] = action['parameters']['source_path']
            action['parameters']['deep'] = True
            self._delete_key(root, lines, action)
            intendation2 = re.search(r'^(\s*)(\S.*)$', lines[dl2])
            StructureChanger.paste_structure(lines, dl2, add, len(intendation2.group(1)) < dc2)        
        return True
        
    def _add_array(self, root, lines, action):
        """Move key transformation"""
        try:
            parent2 = root.get_node_at_path(action['parameters']['addition_path'])
        except:
            return False
        if parent2.implementation != DataNode.Implementation.sequence:
            raise TransformationFileFormatError(
                    "Specified addition path (" + self._get_paths_str(action, 'addition_path') + \
                    ") is not array type node." )
        if parent2.parent is None:
            raise TransformationFileFormatError(
                    "Specified addition path can't be in root node" )
        try:
            parent1 = root.get_node_at_path(action['parameters']['source_path'])
        except:
            return False
        if parent1.implementation != DataNode.Implementation.sequence:
            raise TransformationFileFormatError(
                    "Specified source path (" + self._get_paths_str(action, 'path') + \
                    ") is not array type node." )        
        if parent2 == parent1:
            raise TransformationFileFormatError(
                "Source and addition paths must be different (" + 
                self._get_paths_str(action, 'source_path') + ")")
        sl1, sc1, sl2, sc2 = StructureChanger.value_pos(parent1)
        sl2, sc2 = self._fix_end(lines, sl1, sc1, sl2, sc2 )
        al1, ac1, al2, ac2 = StructureChanger.value_pos(parent2)
        al2, ac2 = self._fix_end(lines, al1, ac1, al2, ac2 )
        if (al1>=sl1 and sl2>=al2) or (sl1>=al1 and al2>=sl2):
            raise TransformationFileFormatError(
                "Source and addition nodes can't contains each other (" + 
                self._get_paths_str(action, 'source_path') + ")")
        
        intendation1 = re.search(r'^(\s*)(\S.*)$', lines[sl1])
        intendation1 = len(intendation1.group(1))
        al1, ac1, al2, ac2 = StructureChanger._add_comments(lines, al1, ac1, al2, ac2)
        add = StructureChanger.copy_structure(lines, al1, ac1, al2, ac2, intendation1)
 
        if al2 < sl1 or (al2 == sl1 and al2 < sc1):
            # source after addition, first delete
            action['parameters']['path'] = action['parameters']['addition_path']
            action['parameters']['deep'] = True
            self._delete_key(root, lines, action)
        self._add_comma(lines, sl2, sc2 )
        StructureChanger.paste_structure(lines, sl2,  add, True)
        if not(al2 < sl1 or (al2 == sl1 and al2 < sc1)):
            action['parameters']['path'] = action['parameters']['addition_path']
            action['parameters']['deep'] = True
            self._delete_key(root, lines, action)
        return True

    def _rename_type(self, root, lines, action):
        """Rename type transformation"""
        try:
            node = root.get_node_at_path(action['parameters']['path'])
        except:
            return False
        old = '!' + action['parameters']['old_name'] 
        new = '!' + action['parameters']['new_name']
 
        StructureChanger.change_tag(lines, node, old,  new)
        return True

    def _fix_end(self, lines, l1, c1, l2, c2 ):
        """If end of node is empty  on next line, end is move to end of preceding line"""
        if c2==0 or lines[l2][:c2].isspace():
            if l1 < l2:
                return l2-1, len(lines[l2-1])
        return l2, c2
        
    def _add_comma(self, lines, l2, c2 ):
        """add comma to the end of array"""
        sep = re.search(r'^(\s*)-(\s+)$', lines[l2])
        if sep is not None:
            return
        if c2 < len(lines[l2]) and not lines[l2][c2:].isspace():
           add =  lines[l2][c2:]
           lines[l2] = lines[l2][:c2]
           lines.insert(l2+1, add)
        lines[l2] += ',' 
       
    def _change_value(self, root, lines, action):
        """Change value to speciffied"""
        try:
            node = root.get_node_at_path(action['parameters']['path'])
        except:
            return False
        if node.implementation != DataNode.Implementation.scalar:
            raise TransformationFileFormatError(
                    "Specified path (" + self._get_paths_str(action, 'path') + ") is not scalar type node." )
        old = '!' + action['parameters']['old_value'] 
        new = '!' + action['parameters']['new_value']
        l1, c1, l2, c2 =  StructureChanger.value_pos(node)
        return StructureChanger.replace(lines, new,  old,  l1, c1, l2, c2 )
        
    def _scale_value(self, root, lines, action):
        """Multiply value by set scale"""
        try:
            node = root.get_node_at_path(action['parameters']['path'])
            scale = float(action['parameters']['scale'])
        except:
            return False
        if node.implementation != DataNode.Implementation.scalar:
            raise TransformationFileFormatError(
                    "Specified path (" + self._get_paths_str(action, 'path') + ") is not scalar type node." )
        try:
            value =  float(node.value)
        except ValueError:
            raise TransformationFileFormatError(
                    "Type of value in specified path (" 
                    + self._get_paths_str(action, 'path') + ") is not numeric." )
        l1, c1, l2, c2 =  StructureChanger.value_pos(node)
        return StructureChanger.replace(lines, str(scale*value),   lines[l1][c1:c2],  l1, c1, l2, c2 )

        
    def _add_key(self, root, lines, action):
        """Add key to abstract record"""
        try:
            node = root.get_node_at_path(action['parameters']['path'])
        except:
            return False
        l1, c1, l2, c2 = StructureChanger.node_pos(node)
        if node.implementation != DataNode.Implementation.mapping:
            raise TransformationFileFormatError(
                "Add path (" + self._get_paths_str(action, 'destination_path') +
                ") must be abstract record")

        key = action['parameters']['key']
        if "value" in action['parameters']:
            value = action['parameters']['value']
        else:
            value = None
        if "type" in action['parameters']:
            tag = action['parameters']['type']
        else:
            tag = None
        if node.parent is None:
            add = StructureChanger.add_key(key, 0, value, tag)
            lines.insert(l1,  add)
        else:
            intendation = re.search(r'^(\s*)(\S.*)$', lines[l1])
            intendation = len(intendation.group(1)) + 2            
            add = StructureChanger.add_key(key, intendation, value, tag)
            lines.insert(l1+1, add)
        return True

class TransformationFileFormatError(Exception):
    """Represents an error in transformation file"""

    def __init__(self, msg):
        super(TransformationFileFormatError, self).__init__(msg)
