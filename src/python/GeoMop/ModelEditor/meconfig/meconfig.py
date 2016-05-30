"""Model Editor static parameters

.. codeauthor:: Pavel Richter <pavel.richter@tul.cz>
.. codeauthor:: Tomas Krizek <tomas.krizek1@tul.cz>
"""

import os
import logging
from copy import deepcopy

import config as cfg
from helpers import (notification_handler, AutocompleteHelper,
                     StructureAnalyzer, shortcuts)
from ist import InfoTextGenerator

from data.import_json import parse_con, fix_tags, rewrite_comments, fix_intendation
from data import export_con
from data.yaml import Loader
from data.yaml import Transformator, TransformationFileFormatError
from data.validation import Validator
from data.format import get_root_input_type_from_json
from data.autoconversion import autoconvert
from geomop_util.logging import LOGGER_PREFIX
from geomop_util import Serializable
from geomop_project import Project, InvalidProject
from util import constants


class _Config:
    """Class for ModelEditor serialization"""

    __serializable__ = Serializable(
        excluded=['observers']
    )

    DEBUG_MODE = False
    """debug mode changes the behaviour"""

    SERIAL_FILE = "ModelEditorData"
    """Serialize class file"""

    COUNT_RECENT_FILES = 5
    """Count of recent files"""

    CONFIG_DIR = os.path.join(cfg.__config_dir__, 'ModelEditor')

    def __init__(self, **kwargs):

        def kw_or_def(key, default=None):
            """Get keyword arg or default value."""
            return kwargs[key] if key in kwargs else default

        from os.path import expanduser
        self.observers = []
        """objects to be notified of changes (currently only used for project)"""
        self.last_data_dir = kw_or_def('last_data_dir', expanduser("~"))
        """directory of the most recently opened data file"""
        self.recent_files = kw_or_def('recent_files', [])
        """a list of recently opened files"""
        self.format_files = kw_or_def('format_files', [])
        """a list of format files"""
        self.display_autocompletion = kw_or_def('display_autocompletion', False)
        """whether to display autocompletion automatically"""
        self.symbol_completion = kw_or_def('symbol_completion', False)
        """whether to automatically complete brackets and array symbols"""
        self.shortcuts = kw_or_def('shortcuts',
                                           deepcopy(shortcuts.DEFAULT_USER_SHORTCUTS))
        """user customizable keyboard shortcuts"""
        self.font = kw_or_def('font', constants.DEFAULT_FONT)
        """text editor font"""
        self._project = kw_or_def('_project')
        self._workspace = kw_or_def('_workspace')

        # initialize project and workspace
        self.workspace = self._workspace
        self.project = self._project

    def update_last_data_dir(self, file_name):
        """Save dir from last used file"""
        project_directory = None
        directory = os.path.dirname(os.path.realpath(file_name))
        if self.workspace is not None and self.project is not None:
            project_dir = os.path.join(self.workspace, self.project)
        if project_directory is None or directory != project_dir:
            self.last_data_dir = directory

    @staticmethod
    def open():
        """Open config from saved file (if exists)."""
        config = cfg.get_config_file(_Config.SERIAL_FILE,
                                     _Config.CONFIG_DIR, cls=_Config)
        if config is None:
            config = _Config()
        return config

    def save(self):
        """Save config data"""
        cfg.save_config_file(self.__class__.SERIAL_FILE, self, self.CONFIG_DIR)

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

    @property
    def data_dir(self):
        """Data directory - either a project dir or the last udes dir."""
        if self.workspace and self.project:
            return os.path.join(self.workspace, self.project)
        else:
            return self.last_data_dir

    @property
    def workspace(self):
        """path to workspace"""
        return self._workspace

    @workspace.setter
    def workspace(self, value):
        if value == '' or value is None:
            self._workspace = None
        if value != self._workspace:
            # close project is workspace is changed
            self.project = None
        self._workspace = value

    @property
    def project(self):
        """name of the project in the workspace"""
        return self._project

    @project.setter
    def project(self, value):
        if value == '' or value is None:
            self._project = None
            Project.current = None
        else:
            self._project = value
            try:
                project = Project.open(self._workspace, self._project)
            except InvalidProject:
                self._project = None
            else:
                Project.current = project
        self.notify_all()

    def notify_all(self):
        """Notify all observers about changes."""
        for observer in self.observers:
            observer.config_changed()


class MEConfig:
    """Static data class"""
    notification_handler = notification_handler
    """error handler for reporting and buffering errors"""
    autocomplete_helper = AutocompleteHelper()
    """helpers for handling autocomplete options in editor"""
    format_files = []
    """Array of format files"""
    transformation_files = []
    """Array of transformation files"""
    curr_format_file = None
    """selected format file"""
    config = _Config.open()
    """Serialized variables"""
    curr_file = None
    """Serialized variables"""
    imported_file_name = None
    """if a file was imported, this is its suggested name"""
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
    resource_dir = os.path.join(os.path.split(
        os.path.dirname(os.path.realpath(__file__)))[0], 'resources')
    """path to a folder containing resources"""
    format_dir = os.path.join(resource_dir, 'format')
    """path to a folder containing IST files"""
    transformation_dir = os.path.join(resource_dir, 'transformation')
    """path to a folder containing transformation files"""
    stylesheet_dir = os.path.join(resource_dir, 'css')
    """path to a folder containing Qt stylesheets"""
    info_text_html_root_dir = os.path.join(resource_dir, 'ist_html')
    """path to a root folder for InfoText"""
    logger = logging.getLogger(LOGGER_PREFIX + constants.CONTEXT_NAME)
    """root context logger"""

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
        for file_name in sorted(listdir(cls.format_dir)):
            if (isfile(join(cls.format_dir, file_name)) and
                    file_name[-5:].lower() == ".json"):
                cls.format_files.append(file_name[:-5])
        # reverse sorting 9 - 0
        cls.format_files = cls.format_files[::-1]

    @classmethod
    def _read_transformation_files(cls):
        """read names of transformation files in format files directory"""
        from os import listdir
        from os.path import isfile, join
        for file_name in listdir(cls.transformation_dir):
            if (isfile(join(cls.transformation_dir, file_name)) and
                    file_name[-5:].lower() == ".json"):
                cls.transformation_files.append(file_name[:-5])

    @classmethod
    def get_curr_format_text(cls):
        """return current format file text"""
        from os.path import join
        file_name = join(cls.format_dir, cls.curr_format_file + ".json")
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
        file_name = join(cls.transformation_dir, file + ".json")
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
        cls.imported_file_name = None

    @classmethod
    def open_file(cls, file_name):
        """
        read file

        return: if file have good format (boolean)
        """
        try:
            with open(file_name, 'r') as file_d:
                cls.document = file_d.read().expandtabs(tabsize=2)
            cls.config.update_last_data_dir(file_name)
            cls.curr_file = file_name
            cls.imported_file_name = None
            cls.config.add_recent_file(file_name, cls.curr_format_file)
            cls.update_format()
            cls.changed = False
            cls.sync_project_for_curr_file()
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
            # find available file name
            base_name = os.path.splitext(os.path.basename(file_name))[0]
            cls.imported_file_name = base_name
            i = 1
            dir_path = cls.config.last_data_dir + os.path.sep
            while os.path.isfile(dir_path + cls.imported_file_name + '.yaml'):
                if i > 999:
                    break
                cls.imported_file_name = "{0}{1:03d}".format(base_name, i)
                i += 1
            cls.imported_file_name = dir_path + cls.imported_file_name + '.yaml'
            cls.curr_file = None
            cls.update()
            cls.document = fix_intendation(cls.document, cls.root)
            cls.update()
            cls.document, need_move_forward = fix_tags(cls.document, cls.root)
            cls.update()
            cls.document = rewrite_comments(con, cls.document, cls.root)
            cls.update()
            data = {'actions': [{'action': 'move-key-forward', 'parameters': {'path': '/system'}},
                                {'action': 'delete-key', 'parameters': {'path': '/system', 'deep': True}}]}
            for path in need_move_forward:
                data['actions'].append({'action': 'move-key-forward', 'parameters': {'path': path}})
            transformator = Transformator(None, data)
            cls.document = transformator.transform(cls.document, cls)
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
    def export_file(cls):
        """Export the current YAML document to CON format.

        :return: text of the exported file
        :rtype: str
        """
        return export_con(cls.root)

    @classmethod
    def open_recent_file(cls, file_name):
        """
        read file from recent files

        return: if file have good format (boolean)
        """
        # If we want to use this code, GUI has to be updated as well.
        # format_file = cls.config.get_format_file(file_name)
        # if format_file is not None:
        #     cls.curr_format_file = format_file
        try:
            with open(file_name, 'r') as file_d:
                cls.document = file_d.read()
            cls.config.update_last_data_dir(file_name)
            cls.curr_file = file_name
            cls.config. add_recent_file(file_name, cls.curr_format_file)
            cls.update_format()
            cls.changed = False
            cls.sync_project_for_curr_file()
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

        # handle parameters
        if (Project.current is not None and
                Project.current.is_abs_path_in_project_dir(cls.curr_file)):
            Project.current.merge_params(cls.validator.params)

        StructureAnalyzer.add_node_info(cls.document, cls.root, cls.notification_handler)
        cls.notifications = cls.notification_handler.notifications

    @classmethod
    def update_format(cls):
        """reread json format file and update node tree"""
        if cls.curr_format_file is None:
            return
        text = cls.get_curr_format_text()
        try:
            cls.root_input_type = get_root_input_type_from_json(text)
        except Exception as e:
            cls._report_error("Can't open format file", e)
        else:
            InfoTextGenerator.init(text)
            cls.autocomplete_helper.create_options(cls.root_input_type)
            cls.update()

    @classmethod
    def save_file(cls):
        """save file"""
        cls.update()
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
        else:
            cls.sync_project_for_curr_file()

    @classmethod
    def save_as(cls, file_name):
        """save file as"""
        cls.update()
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
        else:
            cls.sync_project_for_curr_file()

    @classmethod
    def sync_project_for_curr_file(cls):
        """Write current file and params to a project file."""
        if (Project.current is not None and
                Project.current.is_abs_path_in_project_dir(cls.curr_file)):
            Project.current.merge_params(cls.validator.params)
            params = [param.name for param in cls.validator.params]
            Project.current.add_file(cls.curr_file, params)
            Project.current.save()

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
        if res:
            try:
                cls.document = transformator.transform(cls.document, cls)
                if len(transformator.err)>0:
                    if cls.main_window is not None:
                        cls._report_notify(transformator.err)
                    else:
                        for err in transformator.err:
                            print(err)
            except TransformationFileFormatError as err :
                if cls.main_window is not None:
                    cls._report_error("Transformation format error", err)
                else:
                    raise err
                return
            if transformator.new_version in cls.transformation_files:
                cls.set_current_format_file(transformator.new_version)
            else:
                cls.update()

    @classmethod
    def get_shortcut(cls, name):
        """Locate a keyboard shortcut by its action name.

        :param str name: name of the shortcut
        :return: the assigned shortcut
        :rtype: :py:class:`helpers.keyboard_shortcuts.KeyboardShortcut` or ``None``
        """
        shortcut = None
        if name in shortcuts.SYSTEM_SHORTCUTS:
            shortcut = shortcuts.SYSTEM_SHORTCUTS[name]
        elif name in cls.config.shortcuts:
            shortcut = cls.config.shortcuts[name]
        if shortcut:
            return shortcuts.get_shortcut(shortcut)
        return None

    @classmethod
    def _report_error(cls, mess, err):
        """Report an error with dialog."""
        from geomop_dialogs import GMErrorDialog
        if cls.main_window is not None:
            err_dialog = GMErrorDialog(cls.main_window)
            err_dialog.open_error_dialog(mess, err)
        else:
            raise Exception(mess)
        
    @classmethod
    def _report_notify(cls, errs):
        """Report an error with dialog."""
        from geomop_dialogs import GMErrorDialog
        err_dialog = GMErrorDialog(cls.main_window)
        err_dialog.open_error_report_dialog(errs)
