"""Node analyzer.

.. codeauthor:: Pavel Richter <pavel.richter@tul.cz>
"""

from util import Position
from data import DataNode

from .line_analyzer import LineAnalyzer


class NodeAnalyzer:
    """
    Anayze partial yaml node

    Description: This quick party text analyzing contains node
    specific function.
    """
    
    def __init__(self, lines, node):
        self._lines = lines
        """Array of yaml doc lines"""
        self._node = node
        """analyzed node"""
        
    def get_node_structure_type(self):
        """Get node structure type
        :py:class:`helpers.subyaml.node_analyzer.DataNode.StructureType`"""
        if self._node.implementation == DataNode.Implementation.sequence:
            json_start = self.get_start_inner_json_tag()
            if json_start is not None:
                return DataNode.StructureType.json_array
            return DataNode.StructureType.array
        elif self._node.implementation == DataNode.Implementation.mapping:
                json_start = self.get_start_inner_json_tag()
                if json_start is not None:
                    return DataNode.StructureType.json_dict
                return DataNode.StructureType.dict
        return DataNode.StructureType.scalar
        
    def get_start_inner_json_tag(self):
        """Find start possition of inner json"""
        start_pos = self.get_node_key_end()
        end_pos = self._node.span.start
        if len(self._node.children) > 0:
            end_pos = self._node.children[0].start
        if start_pos is not None and end_pos is not None and \
           start_pos.line <= end_pos.line:
            for i in range(start_pos.line, end_pos.line+1):
                if i>len(self._lines):
                    break
                start=0
                if start_pos.line == i:
                   start = start_pos.column - 1
                end = len(self._lines[i-1])
                if end_pos.line == i:
                    end = end_pos.column - 1
                for char in ["{","[" ]:
                    pos = LineAnalyzer.get_separator(self._lines[i-1],char, start, end)
                    if pos != -1:
                        return Position(i, pos+1)
        return None    
        
    def get_end_inner_json_tag(self):
        """Find start possition of inner json"""
        start_pos = self._node.span.start
        end_pos = self._node.span.end
        if len(self._node.children) > 0:
            start_pos = self._node.children[len(self._node.children)-1].end
        if start_pos is not None and end_pos is not None and \
           start_pos.line <= end_pos.line:
            for i in range(start_pos.line, end_pos.line+1):
                if i>len(self._lines):
                    break
                start=0
                if start_pos.line == i:
                   start = start_pos.column - 1
                end = len(self._lines[i-1])
                if end_pos.line == i:
                    end = end_pos.column - 1
                for char in ["}","]" ]:
                    pos = LineAnalyzer.get_separator(self._lines[i-1],char, start, end)
                    if pos != -1:
                        return Position(i, pos+1)
        return None
    
    def get_node_key_end(self):
        """Get possition after key (or tag, anchor, ... if is not None)"""
        start_pos = self._node.start
        if self._node.key is not None and self._node.key.span is not None:
            start_pos = self._node.key.span.end
            dist = -1
            i = start_pos.line-1
            if start_pos.column > len(self._lines[i]):
                # TODO: can the above condition ever be True? If so, there is probably
                # a mistake in the next line's operator (=+ instead of +=)
                i =+ 1  # this sets the i to 1, instead of incrementing it
                if i > self._node.end.line-1:
                    return start_pos
                line = self._lines[i]
            else:
                line = self._lines[i][start_pos.column-1:]
            while dist == -1: 
                dist = LineAnalyzer.get_after_key_area_end(line)
                if dist is not None:
                    if dist != -1:
                        return Position(i+1, dist)
                    else:
                        i += 1
                        if i > self._node.end.line-1:
                            return self._node.end
                        line = self._lines[i]
                else:
                    return Position(i+1, 1)
        return start_pos
        
    def get_root_node(self):
        """get root node"""
        if self._node is None:
            return None
        node = self._node
        while node.parent is not None:
            node = node.parent
        return node
   
    def get_prev_node(self, line):
        """Get node before line"""
        l=line-1
        root = self.get_root_node()
        if root is None:
            return None
        node = None
        while  node is None:
            if l<0:
                return None
            if not self._lines[l].isspace() and len(self._lines[l])>0: 
                node =  root.get_node_at_position(Position(l+1, len(self._lines[l])-1))
            l -= 1
        return node            
    
    def get_parent_for_unfinished(self, line, index, line_text):
        """return parent node for unfinished node"""
        node = self.get_prev_node(line)
        if node is None:
            return None
        indent = LineAnalyzer.get_indent(line_text)
        if line_text.isspace():
            indent = index
        while True:
            if node.parent is None:
                return node
            prev_indent=LineAnalyzer.get_indent(self._lines[node.start.line-1])            
            if prev_indent < indent and self._node.origin != DataNode.Origin.error:
                if len(self._lines[node.start.line-1]) >= prev_indent+2 and \
                    self._lines[node.start.line-1][prev_indent:prev_indent+2] == "- " and \
                    prev_indent >= indent-2 and node.parent is not None and \
                    node.start.line == node.parent.start.line: 
                        #cursor bellow array character 
                        return node.parent                    
                return node
            node = node.parent 
        return None
