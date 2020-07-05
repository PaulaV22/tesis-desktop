import sys
from PySide import QtCore as QtCore
from PySide import QtGui as QtGui
from project.View.compareView_ui import Ui_mainWindow
import project.Controller.AddController as AC
import project.Controller.DeleteController as DC
import project.Controller.EditController as EC
import project.Controller.CompareController as CC
import project.Controller.AlignController as ALC
import os

class MainWindow(QtGui.QMainWindow, Ui_mainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.setupUi(self)
        self.setWindowTitle("Alelles Searcher")
        self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(self.resourcePath("/View/icon.png"))))
        self.dbFilesTree =QtGui.QTreeView(self)
        self.dbFilesTree.setGeometry(QtCore.QRect(60, 80, 771, 791))
        self.dbFilesTree.setObjectName("dbFilesTree")
        self.widCompare.hide()
        self.widAddDB.hide()
        self.widDeleteDB.hide()
        self.widAlign.hide()
        self.buttonMenuAdd.clicked.connect(self.changeToAdd)
        self.buttonMenuDelete.clicked.connect(self.changeToDelete)
        self.buttonCompare.clicked.connect(self.changeToCompare)
        self.buttonMenuEdit.clicked.connect(self.changeToEdit)
        self.buttonMenuAlign.clicked.connect(self.changeToAlign)
        self.show()
        self.buttonMenuCompare.clicked.connect(self.changeToCompare)
        self.changeToMain()
        self.buttonToCompare.clicked.connect(self.changeToCompare)
        self.buttonToDelete.clicked.connect(self.changeToDelete)
        self.buttonToAdd.clicked.connect(self.changeToAdd)
        self.buttonToEdit.clicked.connect(self.changeToEdit)
        self.buttonToAlign.clicked.connect(self.changeToAlign)

        self.AC = AC.AddController(self)
        self.EC = EC.EditController(self)
        self.DC = DC.DeleteController(self)
        self.CC = CC.CompareController(self)
        self.ALC = ALC.AlignController(self)
        self.AC.configureView()
        self.EC.configureView()
        self.DC.configureView()
        self.CC.configureView()
        self.ALC.configureView()
        self.setDatabases()
        self.dbName= ""
        #self.projectPath = self.AC.getProjectPath()
        self.selectDatabase.currentIndexChanged.connect(self.changeTreeFiles)
        self.setDatabase()

    def resourcePath(self,relative_path):
        """ Get absolute path to resource, works for dev and for PyInstaller """
        base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        output = base_path + relative_path
        return output

    def changeTreeFiles(self):
        self.dbList = self.AC.getDatabases()
        selectedIndex = self.selectDatabase.currentIndex()
        self.dbName = self.dbList[selectedIndex]
        self.AC.setDb(self.dbName)
        self.setDatabase()

    def setDatabases(self):
        self.dbList = self.AC.getDatabases()
        self.selectDeleteDatabase.clear()
        self.selectDatabase.clear()
        for db in self.dbList:
            self.selectDatabase.addItem(db)
            self.selectDeleteDatabase.addItem(db)

    def setDatabase(self):
        selectedDb = ""
        if self.dbName in self.dbList:
            selectedDb = self.resourcePath('\Databases\\' + self.dbName)
        else:
            if len(self.dbList)>0:
                selectedDb = self.resourcePath('\Databases\\'+self.dbList[0])
        if (selectedDb!=""):
            self.filemodel = QtGui.QFileSystemModel()
            self.filemodel.setRootPath(selectedDb)
            self.EC.setDb(self.dbName)
            self.AC.setDb(self.dbName)
            self.DC.setDb(self.dbName)
            self.CC.setDb(self.dbName)

            self.EC.setFileModel(self.filemodel)
            self.dbFilesTree.setModel(self.filemodel)
            indexRoot = self.filemodel.index(self.filemodel.rootPath())
            self.dbFilesTree.setRootIndex(indexRoot)


    def changeToCompare(self):
        self.buttonToAdd.setEnabled(True)
        self.buttonToAdd.show()
        self.buttonToCompare.setEnabled(False)
        self.buttonToCompare.show()
        self.buttonToDelete.setEnabled(True)
        self.buttonToDelete.show()
        self.buttonToEdit.setEnabled(True)
        self.buttonToEdit.show()
        self.buttonToAlign.show()
        self.buttonToAlign.setEnabled(True)
        self.widMainMenu.hide()
        self.widDeleteDB.hide()
        self.widAddDB.hide()
        self.widEditDB.hide()
        self.widAlign.hide()
        self.widCompare.show()

    def changeToEdit(self):
        self.buttonToAdd.setEnabled(True)
        self.buttonToAdd.show()
        self.buttonToCompare.setEnabled(True)
        self.buttonToCompare.show()
        self.buttonToDelete.setEnabled(True)
        self.buttonToDelete.show()
        self.buttonToEdit.setEnabled(False)
        self.buttonToEdit.show()
        self.buttonToAlign.show()
        self.buttonToAlign.setEnabled(True)
        self.widMainMenu.hide()
        self.widDeleteDB.hide()
        self.widAddDB.hide()
        self.widCompare.hide()
        self.widAlign.hide()
        self.widEditDB.show()

    def changeToAdd(self):
        self.widMainMenu.hide()
        self.widDeleteDB.hide()
        self.widCompare.hide()
        self.widAddDB.show()
        self.buttonToAdd.setEnabled(False)
        self.buttonToAdd.show()
        self.buttonToCompare.setEnabled(True)
        self.buttonToCompare.show()
        self.buttonToDelete.setEnabled(True)
        self.buttonToDelete.show()
        self.buttonToEdit.setEnabled(True)
        self.buttonToAlign.show()
        self.buttonToAlign.setEnabled(True)
        self.widEditDB.hide()
        self.widAlign.hide()
        self.buttonToEdit.show()

    def changeToDelete(self):
        self.widMainMenu.hide()
        self.widCompare.hide()
        self.widAddDB.hide()
        self.widDeleteDB.show()
        self.buttonToAdd.setEnabled(True)
        self.buttonToAdd.show()
        self.buttonToCompare.setEnabled(True)
        self.buttonToCompare.show()
        self.buttonToDelete.setEnabled(False)
        self.buttonToDelete.show()
        self.buttonToEdit.setEnabled(True)
        self.buttonToAlign.show()
        self.buttonToAlign.setEnabled(True)
        self.widEditDB.hide()
        self.widAlign.hide()
        self.buttonToEdit.show()

    def changeToMain(self):
        self.widCompare.hide()
        self.widAddDB.hide()
        self.widDeleteDB.hide()
        self.buttonToAdd.setEnabled(True)
        self.buttonToAdd.hide()
        self.buttonToCompare.setEnabled(True)
        self.buttonToCompare.hide()
        self.buttonToDelete.setEnabled(True)
        self.buttonToDelete.hide()
        self.buttonToEdit.setEnabled(True)
        self.buttonToEdit.hide()
        self.buttonToAlign.hide()
        self.buttonToAlign.setEnabled(True)
        self.widEditDB.hide()
        self.widAlign.hide()
        self.widMainMenu.show()


    def changeToAlign(self):
        self.buttonToAdd.setEnabled(True)
        self.buttonToAdd.show()
        self.buttonToCompare.setEnabled(True)
        self.buttonToCompare.show()
        self.buttonToDelete.setEnabled(True)
        self.buttonToDelete.show()
        self.buttonToEdit.setEnabled(True)
        self.buttonToEdit.show()
        self.buttonToAlign.show()
        self.buttonToAlign.setEnabled(False)
        self.widMainMenu.hide()
        self.widDeleteDB.hide()
        self.widAddDB.hide()
        self.widEditDB.hide()
        self.widCompare.hide()
        self.widAlign.show()
