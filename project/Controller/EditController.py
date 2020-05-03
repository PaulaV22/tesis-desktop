import os
from PySide.QtCore import *
from project.model import HaplotypesSearcher as HaplotypesSearcher
from Controller import Controller

class EditController(Controller):

    def __init__(self,window):
        Controller.__init__(self,window)
        self.threadPool = QThreadPool()

    def configureView(self):
        self.window.labelSelectedSeq.show()
        self.window.buttonDeleteSeq.setEnabled(False)
        self.window.newSeqName.setEnabled(True)
        self.window.newSeqName.textChanged.connect(self.checkSeqName)
        self.window.seqContent.setEnabled(True)
        self.window.seqContent.textChanged.connect(self.checkSeqContent)
        self.window.buttonAddSeq.setEnabled(False)
        self.window.progressBar_4.hide()
        self.window.labelEditProcess.hide()
        self.window.dbFilesTree.doubleClicked.connect(self.configure)
        self.window.selectDatabase.currentIndexChanged.connect(self.changeDb)
        self.window.buttonDeleteSeq.clicked.connect(self.deleteSeq)
        self.window.buttonAddSeq.clicked.connect(self.addSeq)
        self.window.progressBar_4.setMinimum(0)
        self.window.progressBar_4.setMaximum(0)
        self.window.progressBar_4.setVisible(False)

    def changeDb(self):
        self.dbList = self.getDatabases()
        selectedIndex = self.window.selectDatabase.currentIndex()
        dbName = self.dbList[selectedIndex]
        if dbName:
            self.path = self.resourcePath('\Databases\\' + dbName)
        else:
            if len(self.dbList) > 0:
                self.path = self.resourcePath('\Databases\\' + self.getDb())
        self.window.labelSelectedDir.setText("Carpeta seleccionada: " + self.path)
        self.checkSeqName()

    def setFileModel(self,filemodel):
        self.filemodel = filemodel

    def configure(self):
        index = self.window.dbFilesTree.currentIndex()
        path = (self.filemodel.filePath(index))
        path2 = path.replace(" ","_")
        self.window.labelEditProcess.hide()
        if os.path.isfile(path2):
            self.toDelete = path2
            self.window.labelSelectedSeq.setText("Secuencia seleccionada: "+path)
            dbPath = self.resourcePath("/" + self.getDb())
            self.window.labelEditProcess.setText(dbPath)
            self.window.buttonDeleteSeq.setEnabled(True)
            self.window.buttonAddSeq.setEnabled(False)
        else:
            self.path = path2
            self.window.labelSelectedDir.setText("Carpeta seleccionada: "+path)
            self.window.labelSelectedSeq.setText("Secuencia seleccionada: ")
            self.window.buttonDeleteSeq.setEnabled(False)
            self.window.newSeqName.setEnabled(True)
            self.window.seqContent.setEnabled(True)
            self.window.buttonAddSeq.setEnabled(True)

    def validSeq(self, seq):
        for base in seq:
            if base not in "ATGCN":
                 return False
        return True

    def checkSeqName(self):
        newSeqPath = self.path+"/"+ self.window.newSeqName.text()+".fa"
        if os.path.exists(newSeqPath):
            self.window.buttonAddSeq.setEnabled(False)
            self.window.labelEditProcess.setText("Ya existe una secuencia con ese nombre en el directorio seleccionado")
            self.window.labelEditProcess.show()
        else:
            self.window.labelEditProcess.setText("")
            self.window.labelEditProcess.hide()
            self.checkSeqContent()

    def checkSeqContent(self):
        self.content = self.window.seqContent.toPlainText()
        self.content =  self.content.replace('\n', '')
        self.content =  self.content.replace('\t', '')
        self.seqName = self.window.newSeqName.text()
        valid = False
        if self.content!="" and self.validSeq(self.content):
            valid = True
        if self.content!="" and self.validSeq(self.content) and self.seqName!= "":
            print("valid seq")
            self.window.buttonAddSeq.setEnabled(True)
            self.window.labelEditProcess.setText("El contenido de la secuencia debe contener las letras A C T G N")

        else:
            if not self.validSeq(self.content):
                self.window.buttonAddSeq.setEnabled(False)
                self.window.labelEditProcess.setText("El contenido de la secuencia debe contener las letras A C T G N")
                self.window.labelEditProcess.show()

    def dbReady(self):
        self.window.progressBar_4.hide()
        self.window.labelEditProcess.setText("Secuencia eliminada correctamente")
        self.window.labelEditProcess.show()
        self.window.buttonDeleteSeq.setEnabled(False)
        self.window.buttonAddSeq.setEnabled(False)


    def deleteSeq(self):
        dbName = self.getDb()
        self.window.progressBar_4.show()
        self.window.progressBar_4.repaint()
        self.window.labelEditProcess.setText("Eliminando secuencia")
        self.window.labelEditProcess.show()
        self.HS = HaplotypesSearcher.HaplotypesSearcher(self.dbName)
        self.HS.signals.deletedSeq.connect(self.dbReady)
        self.HS.setOption("deleteSeq")
        self.HS.setSequenceToDelete(self.toDelete)
        self.threadPool.start(self.HS)

    def newSeqReady(self):
        self.window.progressBar_4.hide()
        self.window.labelEditProcess.setText("Secuencia agregada correctamente")
        self.window.labelEditProcess.show()
        self.window.buttonDeleteSeq.setEnabled(False)
        self.window.buttonAddSeq.setEnabled(False)

    def addSeq(self):
        dbName = self.getDb()
        self.window.progressBar_4.show()
        self.window.progressBar_4.repaint()
        self.window.labelEditProcess.setText("Agregando secuencia")
        self.window.labelEditProcess.show()
        self.HS =HaplotypeSearcher.HaplotypesSearcher(dbName)
        self.HS.setDb(dbName)
        self.HS.signals.addedSeq.connect(self.newSeqReady)
        self.HS.setOption("addSeq")
        self.HS.setAddSeqValues(self.path, self.content, self.seqName)
        self.threadPool.start(self.HS)
