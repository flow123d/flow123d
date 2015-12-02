# -*- coding: utf-8 -*-
"""
GeomMop configuration file parsers

This file contains the parsing functions for configuration files of Flow123d.
Currently supports .con format (as specified by Flow123d manual v1.8.2).
"""

import yaml
import re
from enum import Enum
import json
from collections import OrderedDict
from .data_node import ScalarNode, CompositeNode, NodeOrigin

def _represent_ordereddict(dumper, data):
    value = []

    for item_key, item_value in data.items():
        node_key = dumper.represent_data(item_key)
        node_value = dumper.represent_data(item_value)

        value.append((node_key, node_value))

    return yaml.nodes.MappingNode(u'tag:yaml.org,2002:map', value)
    
def _ordereddict_constructor(loader, node):
    try:
        omap = loader.construct_yaml_omap(node)
        return OrderedDict(*omap)
    except yaml.constructor.ConstructorError:
        return loader.construct_yaml_seq(node)

def parse_con(con):
    """
    Parses a configuration file of Flow123d in .con format with given filename.

    Returns the yaml text structure.
    """
    data = _decode_con(con)
    yaml.add_constructor(u'tag:yaml.org,2002:seq', _ordereddict_constructor)
    yaml.add_representer(OrderedDict, _represent_ordereddict)
    data = yaml.dump(data, default_flow_style=False, indent=2)
    return data


def rewrite_comments(con, yaml, data):
    """Read comments from text (con file) and place it to yaml structure"""
    pattern = re.compile(r"\s?=\s?")
    con = pattern.sub(':', con)
    comments = Comments()
    comments.read_comments_from_con(con, data)
    comments.sorte_by_yaml(yaml, data)
    return comments.write_to_yaml(yaml)


def _decode_con(con):
    """Reads .con format and returns read data in form of dicts and lists."""
    pattern = re.compile(r"//.*")
    con = pattern.sub(r'', con)
    pattern = re.compile(r"/\*.*?\*/", re.DOTALL)
    con = pattern.sub(r'', con)   
    pattern = re.compile(r"^([^\s\[\{,=]+)\s*=\s*")
    con = pattern.sub(r'"\2" : ', con)
    pattern = re.compile(r"([\s\[\{,=])([^\s\[\{,=]+)\s*=\s*")
    con = pattern.sub(r'\1"\2" : ', con)
    return json.loads(con, object_pairs_hook=OrderedDict)

def fix_tags(yaml, root):
    """Replase TYPE and refferences by tags"""
    lines = yaml.splitlines(False)
    add_anchor = {}
    anchor_idx = {}
    del_lines = []
    _traverse_nodes(root, lines, add_anchor, anchor_idx, del_lines)
    new_lines = []
    add_lines = {}
    need_move_forward = []
    for i in add_anchor:
        try:
            anchor = root.get_node_at_path(i)
            col = anchor.span.start.column-1
            line = anchor.span.start.line-1
            if anchor.key.span is not None:
                col = anchor.key.span.end.column-1
                line = anchor.key.span.end.line-1
                if len(lines[line]) > col+1:
                    text = lines[line][col+1:]
                    tag = re.search(r'^(\s+!\S+)', text)
                    if tag is not None:
                        col += len(tag.group(1))
            if len(lines[line]) > col+1:
                if col > 1 and lines[line][col-2:col] == "- ":
                    add_lines[line+1] = col * " " + lines[line][col:]
                    lines[line] = lines[line][:col] + "&anchor" + str(anchor_idx[i])
                else:
                    lines[line] = (lines[line][:col] + " &anchor" + str(anchor_idx[i]) + ' ' +
                                   lines[line][col:])
            else:
                lines[line] += " &anchor" + str(anchor_idx[i])
            # if anchor is in same level as ref and ref is before anchor, rename
            first_ref = None
            for ref_node in add_anchor[i]:
                if ref_node.span.start < anchor.span.start:
                    if first_ref is None or first_ref.span.start > ref_node.span.start:
                        first_ref = ref_node
            if (first_ref is not None and first_ref.parent is not None and
                    first_ref.parent.parent is not None and anchor.parent is not None):
                anchor_path = anchor.parent.absolute_path
                ref_path = first_ref.parent.parent.absolute_path
                if len(anchor_path) <= len(ref_path) and ref_path[:len(anchor_path)] == anchor_path:
                    need_move_forward.append(anchor.absolute_path)
                else:
                    if (first_ref.parent.parent.parent is not None and
                            anchor.parent.parent is not None):
                        anchor_path = anchor.parent.parent.absolute_path
                        ref_path = first_ref.parent.parent.parent.absolute_path
                        if (len(anchor_path) <= len(ref_path) and
                                ref_path[:len(anchor_path)] == anchor_path):
                            need_move_forward.append(anchor.parent.absolute_path)

        except:
            continue
    for i in range(0, len(lines)):
        if i in add_lines:
            new_lines.append(add_lines[i])
        if i not in del_lines:
            new_lines.append(lines[i])
    return "\n".join(new_lines), need_move_forward


def _traverse_nodes(node, lines, add_anchor, anchor_idx, del_lines, i=1):
    """
    Traverse node, recursively call function for children,
    resolve type to tag and resolve references.

    return: array of lines for deleting
    """
    if isinstance(node, CompositeNode):
        for child in node.children:
            if isinstance(child, ScalarNode) and child.key.value == "TYPE":
                del_lines.append(child.key.span.start.line-1)
                lines[node.key.span.start.line-1] += " !" + child.value
            elif isinstance(child, ScalarNode) and child.key.value == "REF":
                del_lines.append(child.key.span.start.line-1)
                if not lines[node.key.span.start.line-1][-1:].isspace():
                    lines[node.key.span.start.line-1] += " "
                lines[node.key.span.start.line-1] += "*anchor" + str(i)
                if child.value not in add_anchor:
                    value = child.value
                    if len(value) > 1 and value[0] == '.':
                        if value[:3] == '../':
                            value = "../" + value
                        else:
                            value = "." + value
                        try:
                            ref_node = child.get_node_at_path(value)
                            value = ref_node.absolute_path
                        except LookupError:
                            pass
                    add_anchor[value] = []
                    anchor_idx[value] = i
                    i += 1
                    add_anchor[value].append(child)
            else:
                if isinstance(child, CompositeNode):
                    i = _traverse_nodes(child, lines, add_anchor, anchor_idx, del_lines, i)
    return i


class Comment:
    """Class for one comment"""
    def __init__(self, path, text, after=False):
        self.path = path
        if len(self.path) > 4 and self.path[-4:] == "TYPE":
            self.path = self.path[:-4]
        if len(self.path) > 3 and self.path[-3:] == "REF":
            self.path = self.path[:-3]
        """comment for key"""
        self.after = after
        """comment after key (else above key)"""
        self.rows = text.splitlines(False)
        """Array of comments"""
        self.pos = None
        """position in associated file"""
        self.inden_line = None
        """line for computation of indentation"""

    def set_yaml_position(self, data):
        """set position in yaml file"""
        node = data.get_node_at_path(self.path)
        if self.after:
            self.pos = node.end
        else:
            self.pos = node.start
        self.inden_line = node.start.line-1


class Comments:
    """Class for comments processed"""
    def __init__(self):
        self.comments = []
        """array of comments"""

    class BlokType(Enum):
        """Type of block"""
        array = 0
        dict = 1

    def read_comments_from_con(self, text, data):
        """reread comments from con file"""
        self.comments = []
        path = ""
        lines = text.splitlines(False)
        col = 0
        line = 0
        comm, col, line = self._read_comments(col, line, lines)
        if comm is not None:
            self.comments.append(Comment('/', comm, False))
        type_, col, line = self._read_start_of_block(col, line, lines)
        if type_ is None:
            return text
        self._traverse_child(col, line, type_, lines, path, data)
        comm, col, line = self._read_comments(col, line, lines, False)
        if comm is not None:
            self.comments.append(Comment('/', comm, True))

    def _traverse_child(self, col, line, type_, lines, path, data):
        """process all data in level and call recursively itself for next children"""
        res, col, line = self._read_end_of_block(col, line, lines, type_)
        if res:
            return col, line
        i = 0
        while line < (len(lines)-1) or col < (len(lines[len(lines)-1])):
            comm, col, line = self._read_comments(col, line, lines)
            if type_ == self.BlokType.dict:
                key, col, line = self._read_key(col, line, lines)
                if key is None:
                    res, col, line = self._read_end_of_block(col, line, lines, type_, False)
                    return col, line
                new_path = self._get_real_path(path, key, data)
                if new_path is None:
                    res, col, line = self._read_end_of_block(col, line, lines, type_, False)
                    return col, line
                if comm is not None:
                    self.comments.append(Comment(new_path, comm, False))
                comm, col, line = self._read_comments(col, line, lines)
            else:
                new_path = self._get_real_path(path, str(i), data)
                i += 1
                if new_path is None:
                    res, col, line = self._read_end_of_block(col, line, lines, type_, False)
                    return col, line
            value, col, line = self._read_value(col, line, lines)
            child_type, col, line = self._read_start_of_block(col, line, lines)
            if comm is not None:
                self.comments.append(Comment(new_path, comm, False))
            if child_type is not None:
                col, line = self._traverse_child(col, line, child_type, lines, new_path, data)
            else:
                value, col, line = self._read_value(col, line, lines)
            comm, col, line = self._read_comments(col, line, lines)
            if comm is not None:
                self.comments.append(Comment(new_path, comm, True))
            res, col, line = self._read_end_of_block(col, line, lines, type_)
            if res:
                comm, col, line = self._read_comments(col, line, lines, True)
                if comm is not None:
                    self.comments.append(Comment(path, comm, True))
                return col, line
            res, col, line = self._read_sep(col, line, lines)
            if res:
                comm, col, line = self._read_comments(col, line, lines, True)
                if comm is not None:
                    self.comments.append(Comment(new_path, comm, True))
            else:
                res, col, line = self._read_end_of_block(col, line, lines, type_, False)
                return col, line

    def sorte_by_yaml(self, yaml, data):
        """sort comments by position in yaml file"""
        for coment in self.comments:
            coment.set_yaml_position(data)
        self.comments = sorted(self.comments, key=lambda comment: comment.pos, reverse=True)
        ok = False
        while not ok:
            ok = True
            for i in range(1, len(self.comments)):
                if self.comments[i].pos.line == self.comments[i-1].pos.line:
                    # first after - reverse order
                    if self.comments[i-1].after and not self.comments[i].after:
                        pom = self.comments[i-1]
                        self.comments[i-1] = self.comments[i]
                        self.comments[i] = pom
                        ok = False
                        break
                    elif self.comments[i].after and \
                           len(self.comments[i-1].path) > len(self.comments[i].path):
                        # first parent - reverse order
                        pom = self.comments[i-1]
                        self.comments[i-1] = self.comments[i]
                        self.comments[i] = pom
                        ok = False
                        break

    def write_to_yaml(self, yaml):
        """return yaml text with comments"""
        lines = yaml.splitlines(False)
        for comment in self.comments:
            intend = re.search(r'^(\s*)(\S.*)$', lines[comment.inden_line])
            if comment.after:
                    if (len(lines[comment.pos.line-1])-1) <= comment.pos.column:
                        if len(lines) <= comment.pos.line:
                            lines.append(intend.group(1) + "# " + comment.rows[0])
                        else:
                            lines.insert(comment.pos.line, intend.group(1) + "# " + comment.rows[0])
                    else:
                        if len(lines) <= comment.pos.line-1:
                            lines.append(intend.group(1) + "# " + comment.rows[0])
                        else:
                            lines.insert(comment.pos.line-1, intend.group(1) + "# " + comment.rows[0])
            else:
                lines.insert(comment.pos.line-1, intend.group(1) + "# " + comment.rows[0])
            for i in range(1, len(comment.rows)):
                if comment.after:
                    lines.insert(comment.pos.line+i, intend.group(1) + "# " + comment.rows[i])
                else:
                    lines.insert(comment.pos.line+i-1, intend.group(1) + "# " + comment.rows[i])
        return "\n".join(lines)

    def _read_key(self, col, line, lines):
        """
        find key in tex

        return key (succes) or None and new line and column
        """
        if line >= len(lines):
            return None, col, line
        txt = lines[line][col:]
        i = 0
        prev_key = None
        icol, iline = self._find_all_interupt(col, line, lines)
        while True:
            key = re.search(r'^(\s*)(\S+)(\s*:)', txt)
            inter = None
            if line + i == iline:
                if i == 0:
                    inter = icol - col
                else:
                    inter = icol
            if key:
                if inter is not None:
                    end_key = len(key.group(1)) + len(key.group(2))
                    if end_key > inter:
                        return None, col, line
                    else:
                        if i == 0:
                            end_key += col
                        end_key += 1
                        if len(lines[line+i]) <= end_key:
                            icol = 0
                            i += 1
                        return self._trim(key.group(2)), end_key, line+i
            else:
                if inter is not None:
                    if prev_key is not None:
                        key = re.search(r'^(\s*):', txt)
                        if key is not None:
                            return self._trim(prev_key), key.span()[1], line+i
                    return None, col, line
            key = re.search(r'^(\s*)$', txt)
            if key is None:
                if prev_key is not None:
                    return None, col, line
                key = re.search(r'^(\s*)(\S*)(\s*)$', txt)
                if key is None:
                    return None, col, line
                else:
                    prev_key = key.group(2)
            i += 1
            if line+i >= len(lines):
                return None, col, line
            txt = lines[line+i]

    def _find_all_interupt(self, col, line, lines):
        """
        find first interapt char

        return its possition
        """
        if line >= len(lines):
            return col, line
        index = None
        ap = None
        text = lines[line][col:]
        i = 0
        for ch in ['[', ']', '{', '}', ',', '//', '/*', '*/', ':']:
            sep = text.find(ch)
            if sep > -1 and (index is None or index > sep):
                index = sep
        for ch in ['"', "'"]:
            sep = text.find(ch)
            if sep > -1 and (ap is None or ap > sep):
                ap = sep
        if index is not None and (ap is None or ap > index):
            if i == 0:
                return index + col, line + i
            else:
                return index, line + i
        if ap is not None:
            chap = text[ap]
            begin = ap+1
            bi = i
            if i == 0:
                begin += col
            if len(lines[i+line]) <= begin:
                i += 1
                begin = 0
                if line+i >= len(lines):
                    return 0, line+i
            text = lines[line+i][begin:]
            while True:
                sep = text.find(chap)
                if sep > -1:
                    new_begin = sep + 1
                    if bi == i:
                        new_begin += begin
                    if len(lines[i+line]) <= new_begin:
                        i += 1
                        new_begin = 0
                        if line+i >= len(lines):
                            return 0, line+i
                    return self._find_all_interupt(new_begin, line+i, lines)
                i += 1
                if line+i >= len(lines):
                    return 0, line+i
                text = lines[line+i]
        i += 1
        if line+i >= len(lines):
            return 0, line+i
        return self._find_all_interupt(0, line+i, lines)

    def _trim(self, text):
        """trim white spaces and apostrophe"""
        if text is None:
            return None
        text = text.strip()
        if len(text) > 1:
            if text[0] == '"' and text[-1] == '"':
                text = text[1:-1]
            if text[0] == "'" and text[-1] == "'":
                text = text[1:-1]
        return text

    def _read_value(self, col, line, lines):
        """
        find key in text

        return value (succes) or line and new line and column
        """
        if line >= len(lines):
            return None, col, line
        txt = lines[line][col:]
        i = 0
        value = None
        icol, iline = self._find_all_interupt(col, line, lines)
        while True:
            inter = None
            if line+i == iline:
                if i == 0:
                    inter = icol - col
                else:
                    inter = icol
            if inter is None:
                if value is None:
                    value = txt.strip()
                else:
                    value += "\n" + txt.strip()
            else:
                if value is None:
                    txt = txt[0:inter]
                else:
                    value += "\n" + txt[0:inter].strip()
                end_col = inter
                if i == 0:
                    end_col += col
                if len(lines[line+i]) <= end_col:
                    end_col = 0
                    i += 1
                return self._trim(value), end_col, line+i
            i += 1
            if line+i >= len(lines):
                return self._trim(value), 0, line+i
            txt = lines[line+i]

    def _read_comments(self, col, line, lines, only_after=False):
        """
        read all comments in text        
        """
        ret  = None
        while True:
            comm, col, line = self._read_comment(col, line, lines, only_after)
            if comm is not None:
                if ret is None:
                    ret = comm
                else:
                    ret += "\n" + comm 
            else:
                return ret, col, line

    def _read_comment(self, col, line, lines, only_after):
        """
        find comment in text

        return comment (success) or None and new line and column
        """
        if only_after and col == 0:
            return None, col, line
        if line >= len(lines):
            return None, col, line
        text = lines[line][col:]
        i = 0
        while True:
            index = None
            for ch in ['//', '/*']:
                sep = text.find(ch)
                if sep > -1 and (index is None or index > sep):
                    index = sep
            if index is not None:
                intend = re.search(r'^(\s*)(\S*)', text)
                if len(intend.group(1)) < index:
                    return None, col, line
                if text[index+1] == '/':
                    if len(text) <= index+2:
                        return None, 0, line+i+1
                    ret = text[index+2:].strip()
                    if len(ret) == 0:
                        return None, 0, line+i+1
                    return ret, 0, line+i+1
                if only_after:
                    return None, col, line
                # /*
                if len(text) <= index+2:
                    text = ""
                else:
                    text = text[index+2:]
                ret = None
                ni = i
                while True:
                    sep = text.find("*/")
                    if sep > -1:
                        if ret is None:
                            ret = text[:sep].rstrip()
                        else:
                            if len(text[:sep].rstrip()) > 0:
                                # last empty line ignore
                                ret += "\n" + text[:sep].rstrip()
                        end = sep+2
                        if i == 0:
                            end += col
                        if i == ni:
                            end += index+2
                        if len(lines[i]) <= end:
                            if len(ret) == 0:
                                return None, 0, line+i+1
                            else:
                                return ret, 0, line+i+1
                        else:
                            if len(ret) == 0:
                                return None, end, line+i
                            else:
                                return ret, end, line+i
                    else:
                        if text.isspace() or len(text) == 0:
                            if ret is None:
                                if i > ni:
                                    # first empty line ignore
                                    ret = "\n"
                            else:
                                ret += "\n"
                        else:
                            if ret is None:
                                ret = text.rstrip()
                            else:
                                ret += "\n" + text.rstrip()
                    i += 1
                    if line+i >= len(lines):
                        return None, 0, line+i
                    text = lines[line+i]
            if len(text) > 0 and not text.isspace():
                return None, col, line
            if only_after:
                return None, col, line
            i += 1
            if line+i >= len(lines):
                return None, col, line
            text = lines[line+i]

    def _read_end_of_block(self, col, line, lines, type_, restricted=True):
        """
        try read end of block, and return its position

        return bool (success) and new line and column
        """
        if line >= len(lines):
            return False, col, line
        if type_ == self.BlokType.dict:
            char = '}'
        else:
            char = ']'
        icol, iline = self._find_all_interupt(col, line, lines)
        if len(lines) <= iline:
            if restricted:
                return False, col, line
            else:
                return False, icol, iline
        if not restricted:
            in_array = 0
            in_dict = 0
            icol, iline = self._find_all_interupt(icol, iline, lines)
            while lines[iline][icol] != char or in_array > 0 or in_dict > 0:
                if len(lines) <= iline:
                    return False, icol, iline
                if lines[iline][icol] == '{':
                    in_dict += 1
                if lines[iline][icol] == '}':
                    in_dict -= 1
                if lines[iline][icol] == '[':
                    in_array += 1
                if lines[iline][icol] == ']':
                    in_array -= 1
                icol += 1
                if len(lines[iline]) <= icol:
                    iline += 1
                    icol = 0
                    if len(lines) <= iline:
                        return True, icol, iline
                icol, iline = self._find_all_interupt(icol, iline, lines)
            return True, icol, iline
        if lines[iline][icol] != char:
            return False, col, line
        text = lines[line][col:]
        i = 0
        while len(text) == 0 or text.isspace():
            i += 1
            if line+i >= len(lines):
                return False, col, line
            text = lines[line+i]
        if iline == (line+i):
            if i == 0:
                text_before = text[:icol-col]
            else:
                text_before = text[:icol]
            if len(text_before) == 0 or text_before.isspace():
                icol += 1
                if len(lines[iline]) <= icol:
                    icol = 0
                    iline += 1
                return True, icol, iline
        return False, col, line

    def _read_sep(self, col, line, lines):
        """
        try read separator, and return its position

        return bool (success) and new line and column
        """
        if line >= len(lines):
            return False, col, line
        icol, iline = self._find_all_interupt(col, line, lines)
        if len(lines) <= iline:
            return False, icol, iline
        if lines[iline][icol] != ',':
            return False, col, line
        text = lines[line][col:]
        i = 0
        while len(text) == 0 or text.isspace():
            i += 1
            if line+i >= len(lines):
                return False, col, line
            text = lines[line+i]
        if iline == (line+i):
            if i == 0:
                text_before = text[:icol-col]
            else:
                text_before = text[:icol]
            if len(text_before) == 0 or text_before.isspace():
                icol += 1
                if len(lines[iline]) <= icol:
                    icol = 0
                    iline += 1
                return True, icol, iline
        return False, col, line

    def _read_start_of_block(self, col, line, lines):
        """
        try read start of block, and return its position

        return type block (success) and new line and column
        """
        if line >= len(lines):
            return None, col, line
        icol, iline = self._find_all_interupt(col, line, lines)
        if len(lines) <= iline:
            return False, icol, iline
        if lines[iline][icol] == '[':
            type_ = self.BlokType.array
        elif lines[iline][icol] == '{':
            type_ = self.BlokType.dict
        else:
            return None, col, line
        text = lines[line][col:]
        i = 0
        while len(text) == 0 or text.isspace():
            i += 1
            if line+i >= len(lines):
                return False, col, line
            text = lines[line+i]
        if iline == (line+i):
            if i == 0:
                text_before = text[:icol-col]
            else:
                text_before = text[:icol]
            if len(text_before) == 0 or text_before.isspace():
                icol += 1
                if len(lines[iline]) <= icol:
                    icol = 0
                    iline += 1
                return type_, icol, iline
        return None, col, line

    def _get_real_path(self, path, key, data):
        """check path existence in data"""        
        if key != "TYPE" and key != "REF":
            path += "/" + key
        node = data
        for key in path.split('/'):
            autoconversion = True
            while autoconversion:
                if not key:
                    # go to next key
                    break
                if not isinstance(node, CompositeNode):
                    return None
                if node.get_child(key) is None:
                    if len(node.children) == 1 and node.children[0].origin != NodeOrigin.structure:
                        node = node.children[0]                        
                    else:
                        return None
                else:
                    node = node.get_child(key)
                    autoconversion = False
        return node.absolute_path
