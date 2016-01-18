"""Icon manimulation library"""
import PyQt5.QtGui as QtGui
import PyQt5.QtCore as QtCore
import os

__icon_dir__ = os.path.join(os.path.dirname(os.path.realpath(__file__)), "icon")

def get_file(name,  size):
    """Get file path to icon from archiv"""
    path = __icon_dir__
    if size == 24:
        path = os.path.join(path, "24x24")
    elif size == 16:
        path = os.path.join(path, "16x16")
    elif size == 32:
        path = os.path.join(path, "32x32")
    elif size == 48:
        path = os.path.join(path, "48x48")
    elif size == 64:
        path = os.path.join(path, "64x64")
    elif size == 128:
        path = os.path.join(path, "128x128")
    return os.path.join(path,name + ".png")

def get_icon(name,  size):
    """Get QIcon object from icon archiv"""
    return QtGui.QIcon(get_file(name,  size))
    
def get_app_icon(name):
    """Get QIcon object from icon archiv"""
    app_icon = QtGui.QIcon()
    app_icon.addFile(get_file(name,  32), QtCore.QSize(32,32))    
    app_icon.addFile(get_file(name,  48), QtCore.QSize(48,48))    
    app_icon.addFile(get_file(name,  64), QtCore.QSize(64,64))
    app_icon.addFile(get_file(name,  64), QtCore.QSize(128,128))
    return app_icon
    
def get_pixmap(name,  size):
    """Get QPixmap object from icon archiv"""
    return QtGui.QPixmap(get_file(name,  size))
