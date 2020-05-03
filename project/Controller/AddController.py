
from Bio import SeqIO
from PySide.QtGui import *
from PySide.QtCore import *
import os
import shutil
from Controller import Controller
from project.model import HaplotypesSearcher as HaplotypesSearcher

class AddController(Controller):

    def __init__(self,window):
        Controller.__init__(self,window)
        self.threadPool = QThreadPool()

    def configureView(self):
        self.fileName=None
        self.window.labelProcess.hide()
        self.window.labelImport.hide()
        self.window.buttonSelect.clicked.connect(self.importFiles)
        self.window.buttonImport.clicked.connect(self.copyFiles)
        self.window.buttonImport.setEnabled(False)
        self.window.inputDbName.textChanged.connect(self.enableButton)
        self.window.progressBar_2.setMinimum(0)
        self.window.progressBar_2.setMaximum(0)
        self.window.progressBar_2.setVisible(False)


    def importFiles(self):
        self.fileName = QFileDialog().getExistingDirectory(None, '',QDir.homePath(), QFileDialog.ShowDirsOnly)
        if self.fileName:
            if (os.path.isdir(self.fileName)):
                self.window.labelImport.setText("Carpeta seleccionada: "+self.fileName)
                self.window.labelImport.show()
                if self.window.inputDbName.text()!="":
                    self.window.buttonImport.setEnabled(True)

    def enableButton(self):
        if self.window.inputDbName.text()!="":
            if not self.fileName is None:
                self.window.buttonImport.setEnabled(True)
                return
        self.window.buttonImport.setEnabled(False)

    def makeDir(self, path):
        if not os.path.exists(path):
            print ("creo el directorio "+ path)
            os.makedirs(path)

    def isFastaFile(self,file, filepath):
        extension = os.path.splitext(file)[1]
        if extension==".fa" or extension==".fasta":
            with open(filepath, "r") as handle:
                fasta = SeqIO.parse(handle, "fasta")
                return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file
        return False

    def copyFiles(self):
        self.dbName = self.window.inputDbName.text()
        dbPath = self.resourcePath("/Databases/"+self.dbName)
        self.makeDir(dbPath)
        self.window.progressBar_2.show()
        self.window.progressBar_2.repaint()
        self.window.labelProcess.setText("Copiando archivos a directorio del programa")
        self.window.labelProcess.show()
        for bases, dirs, files in os.walk(self.fileName):
            for file in files:
                sequenceOrigin = bases + '/' + file
                pathDest = dbPath+bases[len(self.fileName):]
                sequenceDest =  pathDest + '/' + file
                self.makeDir(pathDest)
                if not self.isFastaFile(file,sequenceOrigin):
                    self.window.labelProcess.setText("ERROR: El archivo "+ file+ " no posee el formato Fasta necesario")
                    self.window.labelProcess.show()
                    return
                #print(sequenceOrigin+ "--> "+ sequenceDest)
                if os.path.exists(sequenceOrigin):
                    shutil.copy2(sequenceOrigin, sequenceDest)
                    print("Archivo copiado")
        self.window.labelProcess.setText("Configurando bases de datos en el sistema")
        self.HS = HaplotypesSearcher.HaplotypesSearcher(self.dbName)
        self.HS.setDb(self.dbName)
        self.HS.setNewDb(self.dbName)
        self.HS.signals.database.connect(self.dbReady)
        self.HS.setOption("configuredb")
        self.threadPool.start(self.HS)

    def dbReady(self, status):
        if status=="finished":
            self.window.labelProcess.setText("La base de datos fue agregada correctamente")
            self.window.progressBar_2.hide()
            self.setDatabases()
        else:
            self.window.labelProcess.setText("Error agregando la base de datos. Reintentar.")
            self.window.progressBar_2.hide()
            self.setDatabases()
