"""Structure changer.

.. codeauthor:: Pavel Richter <pavel.richter@tul.cz>
"""
import re
from data import DataNode
from .node_analyzer import NodeAnalyzer

class NodeDescription:
   def __init__(self, type, key=None):
        """init"""
        self.key = key
        """Name of key"""
        self.type = type
        """Type of parent node"""    

class StructureChanger:
    """
    Function for formating and adding yaml text
    Note:
        This class use l1, c1, l2, c2 instad  Position and Span classes.
        When I created it, I wasn't using it. I was thinking about rafactoring,
        but in classes where is StructureChanger used, is using l1, c1, l2, c2
        more comfort. There is a lot various comparation, assignation and 
        other operations of itself variables, that code is shorter and more clearly.
        
        Use  node_pos and value_pos for conversion from span to 
        l1, c1, l2, c2.
    """
    __indent__ = 2
    """File default indentation"""
    __line_length__ = 60
    """File default indentation"""    
    
    def __init__(self):
        """init"""
        pass

    @staticmethod
    def node_pos(node):
        """return position of node in file"""
        return node.start.line-1, node.start.column-1, node.end.line-1, node.end.column-1
    
    @staticmethod
    def value_pos(node):
        """return value position in file of set node """
        return node.span.start.line-1, node.span.start.column-1, \
            node.span.end.line-1, node.span.end.column-1
            
    @staticmethod
    def key_pos(node):
        """return key position in file of set node """
        if node.key is not None and node.key.span is not None:
            return node.key.span.start.line-1, node.key.span.start.column-1, \
                node.key.span.end.line-1, node.key.span.end.column-1
        return node.start.line-1, node.start.column-1, \
            node.start.line-1, node.start.column,
            
    @staticmethod
    def paste_structure(lines, line, add,  after_line, strip_empty=False):
        """
        Paste add structure accoding after_line parameter 
        before or after set possion.
        """        
        for i in range(len(add)-1, -1, -1):
            if not strip_empty or (len(add[i]) > 0 and not add[i].isspace()):
                if after_line:
                    # to end file
                    lines.insert(line+1, add[i])
                else:
                    lines.insert(line, add[i])
                    
    @staticmethod
    def copy_absent_path(root, lines, source_node,  dest_path):
        """
        Return different path (it is not existing path in  dest_path constracted 
        from source_node ) in array of NodeDescription. This structure is suitable 
        for next processing for some function as is paste_absent_path.
        """
        dest = dest_path.split('/')
        if dest[0] == "":
            dest = dest[1:] 
        known_path = ""           
        node = root
        path_exist=True
        while(path_exist):
            path_exist = False
            for child in node.children_keys:
                if child == dest[0]:
                    known_path += "/" + child
                    node = node.get_child(child)
                    dest = dest[1:]
                    path_exist = True
                    break
        node_struct = []
        for i in range(1, len(dest)+1): 
            if source_node is None:
                return [], None            
            if source_node.parent is not None:
                na = NodeAnalyzer(lines, source_node.parent)
                node_type = na. get_node_structure_type()
            else:
                # root element is dictionary
                node_type =  DataNode.StructureType.dict
            if source_node.key is None:
                node_struct.insert(0, NodeDescription(node_type))
            else:
                node_struct.insert(0, NodeDescription(node_type, dest[len(dest)-i]))
            source_node = source_node.parent
        if known_path == "":
            known_path = "/"
        return node_struct, known_path
        
    @staticmethod
    def paste_absent_path(add, node_struct):
        """
        Add node structure to add array and icrease its indentation.
        """
        add_ident = 0
        add_len = len(add)
        prepend_len = 0
        prepend_ident = 0
        add_dash = False        
        # find indetation
        for i in range(0, len(add)):
            place = re.search(r'^(\s*)(\S.*\S)(\s*)$', add[i])
            if place is not None:
                add_ident = len(place.group(1))
                break        
        # prepend (append) structure
        for i in range(0, len(node_struct)):
            if node_struct[i].type == DataNode.StructureType.dict:
                add.insert(prepend_len, (add_ident + prepend_ident)  * " ")
                if  add_dash:
                    add_dash = False
                    add[prepend_len] += "- "
                    prepend_ident += 2
                add[prepend_len] +=  node_struct[i].key +":"
                prepend_len += 1
                prepend_ident += StructureChanger.__indent__
            if node_struct[i].type == DataNode.StructureType.array:
                add_dash = True
            if node_struct[i].type == DataNode.StructureType.json_array:
                add.insert(prepend_len, (add_ident + prepend_ident)  * " ")
                if  add_dash:
                    add_dash = False
                    add[prepend_len] += "- "
                    prepend_ident += 2
                add[prepend_len] +=  "["
                prepend_len += 1
                add.insert(prepend_len+add_len, (add_ident + prepend_ident)  * " " + "]")
                prepend_ident += StructureChanger.__indent__
            if node_struct[i].type == DataNode.StructureType.json_dict:
                add.insert(prepend_len, (add_ident + prepend_ident)  * " ")
                if  add_dash:
                    add_dash = False
                    add[prepend_len] += "- "
                    prepend_ident += 2
                add[prepend_len] +=  "{ " + node_struct[i].key +":"
                prepend_len += 1
                add.insert(prepend_len+add_len, (add_ident + prepend_ident)  * " " + "}")
                prepend_ident += StructureChanger.__indent__
        # add indentation to origin add variable
        for i in range( prepend_len,  prepend_len + add_len):
            if add_dash:
                if len(add[prepend_len]) >  prepend_ident:
                    add[prepend_len] =  (add_ident + prepend_ident)  * " " + "- " + add[i][add_ident:]
                    prepend_ident += 2
                    add_dash = False
            else:
                add[prepend_len] =  prepend_ident * " " +add[prepend_len]        
        return add
    
    @staticmethod
    def change_tag(lines, node, old,  new):
        """change node tag"""        
        l1, c1, l2, c2 = StructureChanger.node_pos(node)
        return StructureChanger._replace(lines, new,  old,  l1, c1, l2, c2 )        
        
    @staticmethod
    def add_tag(lines, node, new):
        """add node tag"""
        l1, c1, l2, c2 = StructureChanger.key_pos(node)
        return StructureChanger._replace(lines, ': ' + new,  ":",  l1, c1, l2, c2 )
        
    
    @staticmethod    
    def replace(lines, new,  old,  l1, c1, l2, c2):
        """replace first occurence of defined value, and return if is len updated"""        
        for i in range(l1, l2+1):
            if i == l1:
                prefix =  lines[i][:c1]
                line =  lines[i][c1:]
            else:
                prefix =  ""                
                line =  lines[i]
            old_str = re.search(re.escape(old), line)
            if old_str is not None:
                if i == l2:
                    if (old_str.end() + len(prefix)):
                        return False
                if old_str.end() == len(line):
                    lines[i] = prefix + line[:old_str.start()] + new
                else:
                    lines[i] = prefix + line[:old_str.start()] + new + line[old_str.end():]
                if i == l2:
                    return len(new) != len(old)
                return False          
        return False       
        
    @staticmethod
    def delete_structure(lines, l1, c1, l2, c2):
        """Delete structure from yaml file (lines array)"""
        if l1 == l2:
            place = re.search(r'^(\s*)(\S.*\S)(\s*)$', lines[l1])
            if place is not None:
                if ((len(place.group(1)) >= c1) and
                        ((len(lines[l1]) - len(place.group(3))) <= c2)):
                    del lines[l1]
                else:
                    lines[l1] = lines[l1][:c1] + lines[l1][c2:]
            else:
                del lines[l1] 
        else:
            # more lines
            place = re.search(r'^(\s*)(\S.*\S)(\s*)$', lines[l2])
            if place is not None:
                if len(place.group(1)) > c2:
                    if (len(lines[l2]) - len(place.group(3))) <= c2:
                        del lines[l2]
                    else:
                        lines[l2] = place.group(1) + lines[l2][c2:]
            else:
                del lines[l2]
            for i in range(l2-1, l1, -1):
                del lines[i]
            place = re.search(r'^(\s*)(\S.*)$', lines[l1])
            if place is not None:
                if len(place.group(1)) >= c1:
                    del lines[l1]
                elif c1 == c2:
                    lines[l1] = lines[l1][:c1] + lines[l1+1][c1:]
                    del lines[l1+1]
                else:
                    lines[l1] = lines[l1][:c1]
            else:
                del lines[l1]

    @staticmethod
    def copy_structure(lines, l1, c1, l2, c2, indent):
        """
        Copy structure from lines to separate array. Structure is
        move by indentation. 
        """
        add = []
        # try add comments
        if l1 == l2:
            add.append(indent*" " + lines[l1][c1:c2])
        else:
            from_line = l1
            if c1 > 0:
                add.append(indent*" " + lines[l1][c1:])
                from_line += 1
            indentation2 = re.search(r'^(\s*)(\S.*)$', lines[l1])
            if indentation2 is None:
                indentation2 = 0
            else:
                indentation2 = len(indentation2.group(1))
            indentation = indent - indentation2
            for i in range(from_line, l2):
                indentation_test = re.search(r'^(\s*)(\S.*)$', lines[i])
                if indentation_test is None:
                    indentation_test = 0
                else:
                    indentation_test = len(indentation_test.group(1))
                if indentation == 0 or indentation_test < -indentation:
                    add.append(lines[i])
                elif indentation < 0:
                    add.append(lines[i][-indentation:])
                else:
                    add.append(indentation*" " + lines[i])
            indentation_test = re.search(r'^(\s*)(\S.*)$', lines[l2])
            if indentation_test is None:
                indentation_test = 0
            else:
                indentation_test = len(indentation_test.group(1))
            if indentation_test < c2:
                if indentation == 0 or indentation_test < -indentation:
                    add.append(lines[l2][:c2])
                elif indentation < 0:
                    add.append(lines[l2][-indentation:c2])
                else:
                    add.append(indentation*" " + lines[l2][:c2])
        return add
        
    @staticmethod
    def skip_tar(lines, l1, c1, l2, c2):
        """
        Return start possition of value from set possition of node. 
        ( Start node structure witout tag, anchor or ref )               
        """
        for char in ["!", r'\*', r'<<:\s*\*', "&"]:
            nl1 = l1
            nc1 = c1
            if lines[nl1][nc1:].isspace() or len(lines[nl1]) <= nc1:
                nl1 += 1
                nc1 = 0
                if nl1 > l2:
                    return l1, c1
                while lines[nl1].isspace() or len(lines[nl1]) == 0:
                    nl1 += 1
                    if nl1 > l2:
                        return l1, c1
            tag = re.search(r'^(\s*' + char + r'\S*)', lines[nl1][nc1:])
            if tag is not None:
                nc1 += len(tag.group(1))
                if len(lines[nl1]) >= nc1:
                    nl1 += 1
                    nc1 = 0
                    if nl1 > l2:
                        return l1, c1
                c1 = nc1
                l1 = nl1
        return l1, c1
        
    @staticmethod
    def _add_comments(lines, l1, c1, l2, c2):
        """Try find comments before and after node and add it to node"""
        nl1 = l1
        nc1 = c1
        nl2 = l2
        nc2 = c2
        # comments before
        inten = re.search(r'^(\s*)(\S+).*$', lines[l1])
        if inten is not None and len(inten.group(1)) >= c1:
            while nl1 > 0:
                comment = re.search(r'^(\s*)#\s*(.*)$', lines[nl1-1])
                if comment is not None:
                    nl1 -= 1
                    nc1 = 0
                else:
                    break
        if nl1 != l1:
            inten = re.search(r'^(\s*)(\S+).*$', lines[nl1])
            nc1 = len(inten.group(1))
        # comments after
        inten = re.search(r'^(.*\S)\s*#\s*\S+.*$', lines[l2])
        if inten is not None and len(inten.group(1)) <= c2:
            nc2 = len(lines[nl2])
        # delete all line comment after
        if c2 == 0:
            comment = re.search(r'^(\s*)#\s*(.*)$', lines[nl2-1])
            if comment is not None:
                nl2 -= 1
                nc2 = len(lines[nl2])
        while nl2 > nl1:
            comment = re.search(r'^(\s*)#\s*(.*)$', lines[nl2])
            if comment is not None:
                nl2 -= 1
                nc2 = len(lines[nl2])
            else:
                break
        return nl1, nc1, nl2, nc2
        
    @staticmethod
    def leave_comments(lines, l1, c1, l2, c2):
        """Try find comments in start and end of node and exclude it from node"""
        nl1 = l1
        nc1 = c1
        nl2 = l2
        nc2 = c2
        # comments on start
        comment = re.search(r'^(\s*)#\s*(.*)$', lines[l1][c1:])
        if comment is not None:
            nl1 += 1
            nc1 = 0
            while nl1 < nl2:
                comment = re.search(r'^(\s*)#\s*(.*)$', lines[nl1])
                if comment is not None:
                    nl1 += 1
                else:
                    break
        while nl2 > nl1:
            comment = re.search(r'^(\s*)#\s*(.*)$', lines[nl2][:nc2])
            if comment is None:
                if nc2 == 0 or lines[nl2][:nc2].isspace():
                    comment = re.search(r'^(\s*)#\s*(.*)$', lines[nl2-1])
                    if comment is None:
                        break
                    else:
                        nl2 -= 1
                        nc2 = 0
                else:
                    nc2 = 0
                ident = re.search(r'^(\s*)\S', lines[nl2])
                if ident is not None:
                    nc2 = len(ident.group(1))
            else:
                nl2 -= 1
                nc2 = 0
                ident = re.search(r'^(\s*)\S', lines[nl2])
                if ident is not None:
                    nc2 = len(ident.group(1))
        return nl1, nc1, nl2, nc2  
      
    @staticmethod
    def add_delete_item_chars(lines, l1, c1, l2, c2):
        """
        If is deleted array node, is needed delete next empty chars and center
        next item as array or delete array char too
        """
        nl1 = l1
        nc1 = c1
        nl2 = l2
        nc2 = c2
        item = re.search(r'^(\s*)-\s(\S.*\S)(\s*)$', lines[l1])
        if item is None:
            item = re.search(r'^(\s*)-\s(\S)(\s*)$', lines[l1])
        if item is not None and \
            (len(item.group(1)) +2) == c1 and \
            (l2>l1 or (len(item.group(1)) + 2 + len(item.group(2)))<=c2):
            # delete array item
            if l1 < l2 and len(lines[l2]) > c2 and not lines[l2][c2:].isspace():
               return l1, c1, l2, c2
            ident = re.search(r'^(\s*)\S', lines[l2+1])
            if ident is not None and len(ident.group(1)) == (len(item.group(1)) + 2):
                #next row as array item
                nl2 = l2 + 1
                nc2 = c1
        return nl1, nc1, nl2, nc2         
            

