"""Module for handling autocomplete in editor.

.. codeauthor:: Tomas Krizek <tomas.krizek1@tul.cz>
"""
from copy import copy
from PyQt5.Qsci import QsciScintilla

from data.format import is_scalar
from geomop_project import Project


class AutocompleteHelper:
    """Creates autocomplete options and manages the autocomplete state."""

    SORTING_ALPHABET = [chr(x) for x in range(48, 57)]  # generate number 0-9
    SORTING_ALPHABET.extend([chr(x) for x in range(97, 123)])  # generate lowercase a-z
    SORTING_ALPHABET.extend(['!', '<', '>', '*', '_', '-'])

    def __init__(self, editor=None):
        """Initialize the class.

        :param YamlEditorWidget editor: instance of editor that can display and hide
           autocomplete
        """
        # TODO: move options to separate class
        self._options = {}
        """the options and their type (anchor, key, ...)"""
        self._anchors = []
        """a list of defined anchors"""
        self.scintilla_options = ''
        """the QScintilla options string encoded in utf-8"""
        self.possible_options = []
        """a list of possible options for the current context"""
        self._visible = False
        """whether autocompletion list is displayed"""
        self.context = None
        """:py:class:`AutocompleteContext` for the current word"""
        self._editor = None
        """editor component that receives signals to show and hide autocompletion"""
        self.pending_check = False
        """is autocompletion pending to be checked for automatic display"""
        self.set_editor(editor)

    @property
    def visible(self):
        """whether autocompletion list is displayed"""
        return self._visible

    @visible.setter
    def visible(self, value):
        """set visibility and emit appropriate signals if state has changed"""
        was_visible = self._visible
        self._visible = value
        if was_visible and not self._visible:
            self._hide()
        elif not was_visible and self._visible:
            self._display()

    def set_editor(self, editor):
        """Set the editor that receives show and hide method calls.

        :param YamlEditorWidget editor: editor that can show and hide autocompletion
        """
        self._editor = editor
        if hasattr(editor, 'SCN_AUTOCCANCELLED'):
            editor.SCN_AUTOCCANCELLED.connect(self._handle_autocompletion_canceled)

    def create_options(self, input_type):
        """Create a list of options based on the input type.

        Each option is identified by a string that should be displayed as QScintilla
        autocomplete option. Set the autocomplete state to hidden.

        :param dict input_type: specification of the input_type
        :return: string of sorted options separated by a space (for QScintilla API)
        :rtype: str
        """
        prev_options = copy(self._options)
        self._options.clear()

        # parameter options
        if Project.current is not None and is_scalar(input_type):
            self._options.update({
                '<' + param.name + '>': 'param'
                for param in Project.current.params
            })

        if input_type['base_type'] == 'Record':  # input type Record
            self._options.update({key: 'key' for key in input_type['keys'] if key != 'TYPE'})
            if 'implemented_abstract_record' in input_type:
                self._options.update({'!' + type_: 'type' for type_ in
                                      input_type['implemented_abstract_record']['implementations']})
        elif input_type['base_type'] == 'Selection':  # input type Selection
            self._options.update({value: 'selection' for value in input_type['values']})
        elif input_type['base_type'] == 'Abstract':  # input typeAbstract
            self._options.update({'!' + type_: 'type' for type_ in
                                  input_type['implementations']})
        elif input_type['base_type'] == 'Bool':  # input tye Bool
            self._options.update({'true': 'bool', 'false': 'bool'})

        # add anchors
        self._options.update({'*' + anchor: 'anchor' for anchor in self._anchors})
        self._prepare_options()

        # if options changed, hide autocomplete
        if sorted(prev_options) != sorted(self._options):
            self.visible = False

        return self.possible_options

    def show_autocompletion(self, context=None):
        """Get autocomplete options for the context.

        If there are some options to be displayed, :py:attr:`visible` is set to True.

        :param AutocompleteContext context: current word and position
        """
        self.refresh_autocompletion(context, create_options=True)
        if len(self.scintilla_options) > 0:
            self.visible = True

    def hide_autocompletion(self):
        """Set the autocompletion state to hidden."""
        self.visible = False

    def refresh_autocompletion(self, context=None, create_options=False):
        """Refresh possible autocomplete option based on context.

        Do not show any options if the autocomplete is hidden. If there are no
        options, set the autocomplete state to hidden.

        :param AutocompleteContext context: current word and position
        :param bool create_options: create options even if autocompletion is not visible
        """
        if not create_options and not self.visible:
            return
        if context is None:
            context = getattr(self._editor, 'autocompletion_context', None)
        self.context = context
        if context.hint is None:
            return
        filter_ = None
        if context is not None:
            filter_ = context.hint
        prev_options = self.scintilla_options
        self._prepare_options(filter_)
        if len(self.scintilla_options) == 0:
            self.visible = False
            return
        if prev_options != self.scintilla_options:
            self._display()
        return

    def get_autocompletion(self, option):
        """Get autocompletion string for the selected option.

        Set the autocompletion state to hidden if an option is selected.

        :param str option: option string returned by QScintilla
        :return: an autocompletion string that replaces the word in current context
        :rtype: str
        :raises: ValueError - if option is not one of possible options
        """
        if option not in self.possible_options:
            raise ValueError("Selected autocompletion option is not available.")
        self.visible = False
        type_ = self._options[option]
        if type_ == 'key':
            return option + ': '
        else:
            return option

    def register_anchor(self, anchor_name):
        """Registers an anchor by its name."""
        self._anchors.append(anchor_name)

    def clear_anchors(self):
        """Clears the anchor list."""
        self._anchors.clear()

    def _prepare_options(self, filter_=None):
        """Sort filtered options and prepare QScintilla string representation.

        :param str filter_: only allow options that starts with this string
        """
        if filter_ is None:
            filter_ = ''
        options = [option for option in self._options.keys() if option.startswith(filter_)]
        self.possible_options = sorted(options, key=self._sorting_key)

        # if there is only one option and it matches the word exactly, do not show options
        if len(self.possible_options) == 1 and self.possible_options[0] == filter_:
            self.possible_options.clear()

        self.scintilla_options = (' '.join(self.possible_options)).encode('utf-8')

    def _sorting_key(self, word):
        """A key for sorting the options."""
        numbers = []
        for letter in word.lower():
            numbers.append(self.SORTING_ALPHABET.index(letter))
        return numbers

    def _display(self):
        """Notify the editor to display autocompletion."""
        send_scintilla = getattr(self._editor, 'SendScintilla', None)
        if callable(send_scintilla):
            send_scintilla(QsciScintilla.SCI_AUTOCSHOW, len(self.context.hint),
                           self.scintilla_options)

    def _hide(self):
        """Notify the editor to hide autocompletion."""
        send_scintilla = getattr(self._editor, 'SendScintilla', None)
        if callable(send_scintilla):
            send_scintilla(QsciScintilla.SCI_AUTOCCANCEL)

    def _handle_autocompletion_canceled(self):
        """Handle when user has cancelled the autocompletion."""
        self._visible = False


class AutocompleteContext:
    """Holds the context information about the word being autocompleted."""

    def __init__(self, word=None, index=None):
        """Initialize the class.

        :param str word: the word being autocompleted
        :param int index: the position of cursor from the beginning of the word
        """
        self.word = word
        """the entire word to be replaced if the autocompletion is triggered"""
        self.index = index
        """the position of cursor from the beginning of the word"""

    @property
    def hint(self):
        """the beginning of the word up to the cursor"""
        if self.word is None or self.index is None:
            return ''
        end = self.index
        return self.word[:end]
