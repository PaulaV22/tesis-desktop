import sys
from project.Controller import *
from project.model import *
from project.View import *
from project.Controller.Main import MainWindow
from PySide import QtCore, QtGui
if __name__ == '__main__':
    __spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"
    try:
        app = QtGui.QApplication(sys.argv)
        mainWin = MainWindow()
        ret = app.exec_()
        sys.exit(ret)
    except Exception as e:
        print(e)