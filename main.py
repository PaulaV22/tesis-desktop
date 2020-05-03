import sys
from project.Controller.Main import MainWindow
from PySide.QtGui import *
from PySide.QtCore import *
from PySide import QtCore, QtGui


if __name__ == '__main__':
    __spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"
    app = QApplication(sys.argv)
    mainWin = MainWindow()
    ret = app.exec_()
    sys.exit(ret)