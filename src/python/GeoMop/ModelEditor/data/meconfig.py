"""Model dialog static parameters"""

import os
import copy

import config as cfg
from helpers import NotificationHandler, AutocompleteHelper
from ist import InfoTextGenerator
from helpers.subyaml import StructureAnalyzer

from .import_json import parse_con, fix_tags, rewrite_comments
from .yaml import Loader
from .yaml import Transformator, TransformationFileFormatError
from .validation import Validator
from .format import get_root_input_type_from_json
from .autoconversion import autoconvert

__author__ = ['Pavel Richter', 'Tomas Krizek']

__resource_dir__ = os.path.join(os.getcwd(), 'resources')
__format_dir__ = os.path.join(__resource_dir__, 'format')
__transformation_dir__ = os.path.join(__resource_dir__, 'transformation')


class _Config:
    """Class for ModelEditor serialization"""

    DEBUG_MODE = False
    """debug mode changes the behaviour"""

    SERIAL_FILE = "ModelEditorData"
    """Serialize class file"""

    COUNT_RECENT_FILES = 5
    """Count of recent files"""

    def __init__(self, readfromconfig=True):
        if readfromconfig:
            data = cfg.get_config_file(self.__class__.SERIAL_FILE)
        else:
            data = None

        if data is not None:
            self.recent_files = copy.deepcopy(data.recent_files)
            self.format_files = copy.deepcopy(data.format_files)
            self.last_data_dir = data.last_data_dir
        else:
            from os.path import expanduser
            self.last_data_dir = expanduser("~")
            self.recent_files = []
            self.format_files = []

    def update_last_data_dir(self, file_name):
        """Save dir from last used file"""
        self.last_data_dir = os.path.dirname(os.path.realpath(file_name))

    def save(self):
        """Save AddPictureWidget data"""
        cfg.save_config_file(self.__class__.SERIAL_FILE, self)

    def add_recent_file(self, file_name, format_file):
        """
        If file is in array, move it top, else add file to top and delete last
        file if is needed. Relevant format files is keep
        """
        # 0 files
        if len(self.recent_files) == 0:
            self.recent_files.append(file_name)
            self.format_files.append(format_file)
            self.save()
            return
        # first file == update file
        if file_name == self.recent_files[0]:
            # format file can be changed
            self.format_files[0] = format_file
            self.save()
            return
        # init for
        last_file = self.recent_files[0]
        last_format = self.format_files[0]
        self.recent_files[0] = file_name
        self.format_files[0] = format_file

        for i in range(1, len(self.recent_files)):
            if file_name == self.recent_files[i]:
                # added file is in list
                self.recent_files[i] = last_file
                self.format_files[i] = last_format
                self.save()
                return
            last_file_pom = self.recent_files[i]
            last_format_pom = self.format_files[i]
            self.recent_files[i] = last_file
            self.format_files[i] = last_format
            last_file = last_file_pom
            last_format = last_format_pom
            # recent files is max+1, but first is not displayed
            if self.__class__.COUNT_RECENT_FILES < i+1:
                self.save()
                return
        # add last file
        self.recent_files.append(last_file)
        self.format_files.append(last_format)
        self.save()

    def get_format_file(self, file_name):
        """get format file that is in same position as file"""
        for i in range(0, len(self.recent_files)):
            if self.recent_files[i] == file_name:
                return self.format_files[i]
        return None


class MEConfig:
    """Static data class"""
    notification_handler = NotificationHandler()
    """error handler for reporting and buffering errors"""
    autocomplete_helper = AutocompleteHelper()
    """helpers for handling autocomplete options in editor"""
    format_files = []
    """Array of format files"""
    transformation_files = []
    """Array of transformation files"""
    curr_format_file = None
    """selected format file"""
    config = _Config()
    """Serialized variables"""
    curr_file = None
    """Serialized variables"""
    root = None
    """root DataNode structure"""
    document = ""
    """text set by editor after significant changing"""
    main_window = None
    """parent of dialogs"""
    notifications = []
    """array of notifications"""
    changed = False
    """is file changed"""
    loader = Loader(notification_handler)
    """loader of YAML files"""
    validator = Validator(notification_handler)
    """data validator"""
    root_input_type = None
    """input type of the whole tree, parsed from format"""

    def __init__(self):
        pass

    @classmethod
    def init(cls, main_window):
        """Init class wit static method"""
        cls._read_format_files()
        cls._read_transformation_files()
        cls.main_window = main_window
        if len(cls.config.format_files) > 0:
            cls.curr_format_file = cls.config.format_files[0]
        else:
            if len(cls.format_files) > 0:
                cls.curr_format_file = cls.format_files[0]
        cls.update_format()

    @classmethod
    def _read_format_files(cls):
        """read names of format files in format files directory"""
        from os import listdir
        from os.path import isfile, join
        for file_name in sorted(listdir(__format_dir__)):
            if (isfile(join(__format_dir__, file_name)) and
                    file_name[-5:].lower() == ".json"):
                cls.format_files.append(file_name[:-5])

    @classmethod
    def _read_transformation_files(cls):
        """read names of transformation files in format files directory"""
        from os import listdir
        from os.path import isfile, join
        for file_name in listdir(__transformation_dir__):
            if (isfile(join(__transformation_dir__, file_name)) and
                    file_name[-5:].lower() == ".json"):
                cls.transformation_files.append(file_name[:-5])

    @classmethod
    def get_curr_format_text(cls):
        """return current format file text"""
        from os.path import join
        file_name = join(__format_dir__, cls.curr_format_file + ".json")
        try:
            with open(file_name, 'r') as file_d:
                return file_d.read()
        except (RuntimeError, IOError) as err:
            if cls.main_window is not None:
                cls._report_error(
                    "Can't open format file '" + cls.curr_format_file + "'", err)
            else:
                raise err
        return None

    @classmethod
    def get_transformation_text(cls, file):
        """return transformation file text"""
        from os.path import join
        file_name = join(__transformation_dir__, file + ".json")
        try:
            with open(file_name, 'r') as file_d:
                return file_d.read()
        except (RuntimeError, IOError) as err:
            if cls.main_window is not None:
                cls._report_error(
                    "Can't open transformation file '" + file + "'", err)
            else:
                raise err
        return None

    @classmethod
    def get_data_node(cls, position):
        """
        Returns DataNode at given `class::Position` position.
        """
        # empty file with comment
        if cls.root is None:
            return None

        return cls.root.get_node_at_position(position)

    @classmethod
    def set_current_format_file(cls, file_name):
        """
        set current format file
        """
        if file_name not in cls.format_files:
            return
        cls.curr_format_file = file_name
        cls.update_format()

    @classmethod
    def new_file(cls):
        """
       empty file
        """
        cls.document = ""
        cls.update_format()
        cls.changed = False
        cls.curr_file = None

    @classmethod
    def open_file(cls, file_name):
        """
        read file

        return: if file have good format (boolean)
        """
        try:
            with open(file_name, 'r') as file_d:
                cls.document = file_d.read()
            cls.config.update_last_data_dir(file_name)
            cls.curr_file = file_name
            cls.config.add_recent_file(file_name, cls.curr_format_file)
            cls.update_format()
            cls.changed = False
            return True
        except (RuntimeError, IOError) as err:
            if cls.main_window is not None:
                cls._report_error("Can't open file", err)
            else:
                raise err
        return False

    @classmethod
    def import_file(cls, file_name):
        """
        read con file and transform it to yaml structure

        return: if file have good format (boolean)
        """
        try:
            with open(file_name, 'r') as file_d:
                con = file_d.read()
            cls.document = parse_con(con)
            cls.curr_file = None
            cls.update()
            cls.document, need_move_forward = fix_tags(cls.document, cls.root)
            cls.update()
            cls.document = rewrite_comments(con, cls.document, cls.root)
            cls.update()
            data = {'actions': [{'action': 'move-key-forward', 'parameters': {'path': '/system'}},
                                {'action': 'delete-key', 'parameters': {'path': '/system'}}]}
            for path in need_move_forward:
                data['actions'].append({'action': 'move-key-forward', 'parameters': {'path': path}})
            transformator = Transformator(None, data)
            cls.document = transformator.transform(cls.document)
            cls.update_format()
            cls.changed = True
            return True
        except (RuntimeError, IOError) as err:
            if cls.main_window is not None:
                cls._report_error("Can't open import file", err)
            else:
                raise err
        except Exception as err:
            if cls.main_window is not None:
                cls._report_error("Can't import file from con format", err)
            else:
                raise err
        return False

    @classmethod
    def open_recent_file(cls, file_name):
        """
        read file from recent files

        return: if file have good format (boolean)
        """
        format_file = cls.config.get_format_file(file_name)
        if format_file is not None:
            cls.curr_format_file = format_file
        try:
            with open(file_name, 'r') as file_d:
                cls.document = file_d.read()
            cls.config.update_last_data_dir(file_name)
            cls.curr_file = file_name
            cls.config. add_recent_file(file_name, cls.curr_format_file)
            cls.update_format()
            cls.changed = False
            return True
        except (RuntimeError, IOError) as err:
            if cls.main_window is not None:
                cls._report_error("Can't open file", err)
            else:
                raise err
        return False

    @classmethod
    def update(cls):
        """reread yaml text and update node tree"""
        cls.notification_handler.clear()
        cls.root = cls.loader.load(cls.document)
        cls.autocomplete_helper.clear_anchors()
        for anchor in cls.loader.anchors:
            cls.autocomplete_helper.register_anchor(anchor)
        cls.notifications = cls.notification_handler.notifications
        if cls.root_input_type is None or cls.root is None:
            return
        cls.root = autoconvert(cls.root, cls.root_input_type)
        cls.validator.validate(cls.root, cls.root_input_type)
        StructureAnalyzer.add_node_info(cls.document, cls.root, cls.notification_handler)
        cls.notifications = cls.notification_handler.notifications

    @classmethod
    def update_format(cls):
        """reread json format file and update node tree"""
        if cls.curr_format_file is None:
            return
        text = cls.get_curr_format_text()
        cls.root_input_type = get_root_input_type_from_json(text)
        InfoTextGenerator.init(text)
        cls.autocomplete_helper.create_options(cls.root_input_type)
        cls.update()

    @classmethod
    def save_file(cls):
        """save file"""
        try:
            file_d = open(cls.curr_file, 'w')
            file_d.write(cls.document)
            file_d.close()
            # format is save to recent files up to save file
            cls.config.format_files[0] = cls.curr_format_file
            cls.changed = False
        except (RuntimeError, IOError) as err:
            if cls.main_window is not None:
                cls._report_error("Can't save file", err)
            else:
                raise err

    @classmethod
    def save_as(cls, file_name):
        """save file as"""
        try:
            file_d = open(file_name, 'w')
            file_d.write(cls.document)
            file_d.close()
            cls.config.update_last_data_dir(file_name)
            cls.curr_file = file_name
            cls.config.add_recent_file(file_name, cls.curr_format_file)
            cls.changed = False
        except (RuntimeError, IOError) as err:
            if cls.main_window is not None:
                cls._report_error("Can't save file", err)
            else:
                raise err

    @classmethod
    def update_yaml_file(cls, new_yaml_text):
        """update new editor text"""
        if new_yaml_text != cls.document:
            cls.document = new_yaml_text
            cls.changed = True
            return True
        return False

    @classmethod
    def transform(cls, file):
        """Run transformation according rules in set file"""
        cls.update()
        text = cls.get_transformation_text(file)
        try:
            transformator = Transformator(text)
        except (ValueError, TransformationFileFormatError) as err:
            if cls.main_window is not None:
                cls._report_error("Can't decode transformation file", err)
            else:
                raise err
            return
        dialog = None
        res = True
        if cls.main_window is not None:
            import PyQt5.QtWidgets as QtWidgets
            from ui.dialogs import TranformationDetailDlg
            
            dialog = TranformationDetailDlg(transformator.name,
                                            transformator.description,
                                            transformator.old_version,
                                            cls.curr_format_file,
                                            transformator.new_version,
                                            transformator.new_version in cls.transformation_files,
                                            cls.main_window)
            res = QtWidgets.QDialog.Accepted == dialog.exec_()
        if res :
            try:
                cls.document = transformator.transform(cls.document)
            except TransformationFileFormatError as err:
                if cls.main_window is not None:
                    cls._report_error("Can't decode transformation file", err)
                else:
                    raise err
                return
            if transformator.new_version in cls.transformation_files:
                cls.set_current_format_file(transformator.new_version)
            else:
                cls.update()
                
    @classmethod
    def _report_error(cls, mess,  err):
        from geomop_dialogs import GMErrorDialog
        err_dialog = GMErrorDialog(cls.main_window)
        err_dialog.open_error_dialog(mess, err)
