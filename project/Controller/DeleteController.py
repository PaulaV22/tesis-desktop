from PySide.QtCore import *
from Controller import Controller
from project.model import HaplotypesSearcher as HaplotypeSearcher

class DeleteController(Controller):

    def __init__(self,window):
        Controller.__init__(self,window)
        self.threadPool = QThreadPool()

    def configureView(self):
        self.selectedDbName=None
        self.window.labelProcessDelete.hide()
        self.window.askingWidgetDelete.hide()
        self.window.buttonDelete.setEnabled(True)
        self.window.buttonDelete.clicked.connect(self.showWidget)
        self.window.progressBar_3.setMinimum(0)
        self.window.progressBar_3.setMaximum(0)
        self.window.progressBar_3.setVisible(False)
        self.window.buttonBox.accepted.connect(self.deleteDb)
        self.window.buttonBox.rejected.connect(self.hideWidget)
        self.setDatabases()


    def setDatabasesOptions(self):
        self.dbList = self.getDatabases()
        for db in self.dbList:
            self.window.selectDeleteDatabase.addItem(db)

    def removeItem(self):
        self.dbList = self.getDatabases()
        index = self.window.selectDeleteDatabase.findText(self.selectedDbName)
        self.window.selectDeleteDatabase.removeItem(index)
        self.window.selectDatabase.removeItem(index)


    def showWidget(self):
        self.window.labelProcessDelete.hide()
        selectedIndex = self.window.selectDeleteDatabase.currentIndex()
        self.selectedDbName = self.dbList[selectedIndex]
        if (self.selectedDbName):
            self.window.askingWidgetDelete.show()

    def hideWidget(self):
        self.window.askingWidgetDelete.hide()

    def deleteDb(self):
        self.window.labelProcessDelete.setText("Eliminando la base de datos")
        self.window.labelProcessDelete.show()
        self.window.progressBar_3.show()
        self.HS = HaplotypeSearcher.HaplotypesSearcher(self.selectedDbName)
        self.HS.signals.deleted.connect(self.deletedDb)
        self.HS.setOption("deletedb")
        self.threadPool.start(self.HS)

    def deletedDb(self):
        self.hideWidget()
        self.window.labelProcessDelete.setText("La base de datos ha sido eliminada correctamente ")
        self.window.progressBar_3.hide()
        self.setDatabases()

