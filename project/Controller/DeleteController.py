from PySide.QtCore import *
from Controller import Controller
from project.model import HaplotypesSearcher as HaplotypeSearcher
from project.model import DbAdmin as DbAdmin

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
        self.dbAdmin = DbAdmin.DbAdmin(self.selectedDbName)
        self.dbAdmin.signals.deleted.connect(self.deletedDb)
        self.dbAdmin.setOption("deletedb")
        self.threadPool.start(self.dbAdmin)

    def deletedDb(self):
        self.hideWidget()
        self.window.labelProcessDelete.setText("La base de datos ha sido eliminada correctamente ")
        self.window.progressBar_3.hide()
        self.setDatabases()

