"""Module for exporting the data structure to .con format.

.. codeauthor:: Tomas Krizek <tomas.krizek1@tul.cz>
"""

INDENTATION = '  '


class Exporter:
    """Exporter from data structure to con files."""

    def __init__(self):
        """Initialize the class."""
        self.lines = ['']
        self.indent_level = 0

    def export_con(self, root):
        """Create .con text from a root data node.

        :param DataNode root: the root of the data structure
        :return: text representation of the structure in .con format
        :rtype: str
        """
        self.lines = ['']
        self.indent_level = 0

        self._create_node(root)

        return '\n'.join(self.lines)

    def _print_line(self, text):
        """Append the text as indented line to the buffer.

        :param str text: a line of text without the EOL symbol
        """
        self.lines.append(self.indent_level * INDENTATION + text)

    def _print(self, text):
        """Append the text to the last line."""
        self.lines[-1] = self.lines[-1] + text

    def _print_new_line(self, indent_change=0):
        """Append new line with the appropriate indentation.

        :param int indent_change: +1, 0 or -1 to increase, keep or decrease indentation
        """
        self.indent_level += indent_change
        self.lines.append(self.indent_level * INDENTATION)

    def _create_mapping_node(self, node):
        """Create a mapping node."""
        self._print('{')
        self._print_new_line(1)

        # check for type
        if node.type is not None:
            self._print('TYPE = "{type}",'.format(type=node.type.value))
            self._print_new_line()

        # print all keys
        for child in node.children:
            self._print(child.key.value + ' = ')
            self._create_node(child)
            self._print(',')
            self._print_new_line()

        self.lines.pop()  # remove last (extra) line
        self.lines[-1] = self.lines[-1][:-1]  # remove , from end of line

        self._print_new_line(-1)
        self._print('}')

    def _create_node(self, node):
        """Create a node based on its type.

        :param DataNode node: node to be create in text
        """
        if node.ref is not None:
            path = node.ref.absolute_path
            self._create_reference(path)
        else:
            if node.implementation == node.Implementation.mapping:
                self._create_mapping_node(node)
            elif node.implementation == node.Implementation.scalar:
                self._create_scalar_node(node)
            elif node.implementation == node.Implementation.sequence:
                self._create_sequence_node(node)

    def _create_scalar_node(self, node):
        """Create a text representation of scalar node.

        :param DataNode node: node
        """
        if isinstance(node.value, bool):
            self._print('true' if node.value else 'false')
        elif isinstance(node.value, int):
            self._print(str(node.value))
        elif isinstance(node.value, float):
            self._print(str(node.value))
        else:
            self._print('"' + node.value + '"')

    def _create_sequence_node(self, node):
        """Create a text representation of sequence node.

        :param DataNode node: node
        """
        self._print('[')
        self._print_new_line(1)

        # print all keys
        for child in node.children:
            self._create_node(child)
            self._print(',')
            self._print_new_line()

        self.lines.pop()  # remove last (extra) line
        self.lines[-1] = self.lines[-1][:-1]  # remove , from end of line

        self._print_new_line(-1)
        self._print(']')

    def _create_reference(self, path):
        """Create a reference node with the given absolute path."""
        self._print('{')
        self._print_new_line(1)
        self._print('REF = "{ref}"'.format(ref=path))
        self._print_new_line(-1)
        self._print('}')


_exporter = Exporter()
export_con = _exporter.export_con
