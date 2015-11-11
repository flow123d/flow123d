"""Module for handling autocomplete in editor."""

__author__ = 'Tomas Krizek'


class AutocompleteHelper:
    """Helper class for creating and managing autocomplete options in editor."""

    def __init__(self):
        """Initializes the class."""
        self._options = {}
        self._anchors = []
        self.scintilla_options = ''
        """the QScintilla options string encoded in utf-8"""

        # define sorting alphabet
        self.sorting_alphabet = list(map(chr, range(48, 57)))  # generate number 0-9
        self.sorting_alphabet.extend(list(map(chr, range(97, 123))))  # generate lowercase alphabet
        self.sorting_alphabet.extend(['!', '*', '_', '-'])

    def create_options(self, input_type):
        """
        Creates a list of options based on the input type. Each option is identified by a string.

        Returns a list of string (option identifiers) that should be displayed as QScintilla
        autocomplete options.
        """
        self._options.clear()

        if input_type['base_type'] == 'Record':  # input type Record
            self._options.update({key: 'key' for key in input_type['keys'] if key != 'TYPE'})
            if 'implemented_abstract_record' in input_type:
                self._options.update({'!' + type_: 'type' for type_ in
                                     input_type['implemented_abstract_record']['implementations']})
        elif input_type['base_type'] == 'Selection':  # input type Selection
            self._options.update({value: 'selection' for value in input_type['values']})
        elif input_type['base_type'] == 'AbstractRecord':  # input typeAbstractRecord
            self._options.update({'!' + type_: 'type' for type_ in
                                  input_type['implementations']})
        elif input_type['base_type'] == 'Bool':  # input tye Bool
            self._options.update({'true': 'bool', 'false': 'bool'})

        # add anchors
        self._options.update({'*' + anchor: 'anchor' for anchor in self._anchors})

        # sort and create scintilla options string
        sorted_options = sorted(list(self._options.keys()), key=self._sorting_key)
        self.scintilla_options = (' '.join(sorted_options)).encode('utf-8')

        return sorted_options

    def select_option(self, option_string):
        """
        Selects an option based on the option_string returned from QScintilla API.

        Returns a string that should be inserted into the editor. This string can differ from
        `option_string` for a better user experience.(For example, a record key would be identified
        by its name, but the returned string would also contain ': ' at the end.
        """
        if option_string not in self._options:
            return ''
        type_ = self._options[option_string]
        if type_ == 'key':
            return option_string + ': '
        else:
            return option_string

    def register_anchor(self, anchor_name):
        """Registers an anchor by its name."""
        self._anchors.append(anchor_name)

    def clear_anchors(self):
        """Clears the anchor list."""
        self._anchors.clear()

    def _sorting_key(self, word):
        """A key for sorting the options."""
        numbers = []
        for letter in word.lower():
            numbers.append(self.sorting_alphabet.index(letter))
        return numbers
