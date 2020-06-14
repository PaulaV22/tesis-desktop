# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'compareView.ui'
#
# Created: Mon Aug 19 23:13:09 2019
#      by: pyside-uic 0.2.15 running on PySide 1.2.4
#
# WARNING! All changes made in this file will be lost!

from PySide import QtCore, QtGui

class Ui_mainWindow(object):
    def setupUi(self, mainWindow):
        mainWindow.setObjectName("mainWindow")
        mainWindow.resize(842, 550)
        self.centralwidget = QtGui.QWidget(mainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.labelDatabase = QtGui.QLabel(self.centralwidget)
        self.labelDatabase.setGeometry(QtCore.QRect(30, 20, 141, 21))
        self.labelDatabase.setObjectName("labelDatabase")
        self.selectDatabase = QtGui.QComboBox(self.centralwidget)
        self.selectDatabase.setGeometry(QtCore.QRect(150, 20, 181, 27))
        self.selectDatabase.setObjectName("selectDatabase")
        self.labelFiles = QtGui.QLabel(self.centralwidget)
        self.labelFiles.setGeometry(QtCore.QRect(30, 60, 70, 21))
        self.labelFiles.setObjectName("labelFiles")
        self.listFiles = QtGui.QListView(self.centralwidget)
        self.listFiles.setGeometry(QtCore.QRect(30, 100, 301, 371))
        self.listFiles.setObjectName("listFiles")
        self.plainTextEdit = QtGui.QPlainTextEdit(self.centralwidget)
        self.plainTextEdit.setGeometry(QtCore.QRect(440, 60, 381, 161))
        self.plainTextEdit.setObjectName("plainTextEdit")
        self.label = QtGui.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(480, 20, 251, 21))
        self.label.setObjectName("label")
        self.buttonCompare = QtGui.QPushButton(self.centralwidget)
        self.buttonCompare.setGeometry(QtCore.QRect(570, 230, 112, 34))
        self.buttonCompare.setObjectName("buttonCompare")
        self.tableResults = QtGui.QTableView(self.centralwidget)
        self.tableResults.setGeometry(QtCore.QRect(440, 280, 381, 192))
        self.tableResults.setObjectName("tableResults")
        self.buttonDownload = QtGui.QToolButton(self.centralwidget)
        self.buttonDownload.setGeometry(QtCore.QRect(390, 440, 32, 27))
        self.buttonDownload.setObjectName("buttonDownload")
        mainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(mainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 842, 31))
        self.menubar.setMouseTracking(True)
        self.menubar.setObjectName("menubar")
        mainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(mainWindow)
        self.statusbar.setObjectName("statusbar")
        mainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(mainWindow)
        QtCore.QMetaObject.connectSlotsByName(mainWindow)

    def retranslateUi(self, mainWindow):
        mainWindow.setWindowTitle(QtGui.QApplication.translate("mainWindow", "MainWindow", None, QtGui.QApplication.UnicodeUTF8))
        self.labelDatabase.setText(QtGui.QApplication.translate("mainWindow", "Base de datos", None, QtGui.QApplication.UnicodeUTF8))
        self.labelFiles.setText(QtGui.QApplication.translate("mainWindow", "Archivos", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("mainWindow", "Ingrese la secuencia a comparar", None, QtGui.QApplication.UnicodeUTF8))
        self.buttonCompare.setText(QtGui.QApplication.translate("mainWindow", "Comparar", None, QtGui.QApplication.UnicodeUTF8))
        self.buttonDownload.setText(QtGui.QApplication.translate("mainWindow", "...", None, QtGui.QApplication.UnicodeUTF8))

