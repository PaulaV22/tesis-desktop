import SimpleDbCreator as SC
import SimpleBlast as S
from project.model import AmbiguousDbCreator as AC, GlobalBlast as GC
import ResultsAnalizer as RA
import os
import shutil
import json
from Signal import HaploSignal
import sys
from PySide.QtGui import *
from PySide.QtCore import *

class DbAdmin(QRunnable):

    def __init__(self, dbName=None):
        QRunnable.__init__(self)
        if (dbName):
            self.db = dbName
        else:
            self.db = "BoLa"
        self.projectPath = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.setDb(self.db)
        self.signals = HaploSignal()
        self.option = "compare"
        self.newDb = None

    def resourcePath(self,relative_path):
        """ Get absolute path to resource, works for dev and for PyInstaller """
        base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        output = base_path + relative_path
        return output

    def setOption(self,option):
        self.option=option

    def setNewDb(self,newDb):
        self.newDb=newDb

    def run(self):
    #def searchHaplotypes(self):
        if self.option=="configuredb":
            try:
                self.configureDb(self.newDb)
                self.signals.database.emit("finished")
                self.signals.updatedatabases.emit()
            except:
                self.signals.database.emit("error")
        if self.option=="deletedb":
            self.deleteDb(self.db)
            self.signals.deleted.emit()
        if self.option=="deleteSeq":
            self.deleteSeq(self.db, self.sequenceToDelete)
            self.signals.deletedSeq.emit()
        if self.option=="addSeq":
            self.addSeq(self.path,self.db,self.newSeqName, self.newSeqContent)
            self.signals.addedSeq.emit()

    def configureDb(self, db):
        ####crear la bd con los archivos originales de BoLa####
        ready = False
        self.simpleDbCreator.makeDb()
        ####alinear todas las secuencias de BoLa entre si generando un archivo de salida por cada alineacion (n x n)####
        while not ready:
            try:
                self.globalBlast.align("/Databases/" + db)
                ready = True
            except:
                self.simpleDbCreator.makeDb()
                ready = False
        ####armar la base de datos con las posibles combinaciones (Nuevadb)####
        self.ambiguousDbCreator.makeDb()

    def deleteDb(self,db,total=True):
        Database = self.resourcePath('/Databases/' + db)
        BlastDb = self.resourcePath('/BlastDb/' + db)
        BlastResult = self.resourcePath('/BlastResult/' + db)
        DbAmbigua = self.resourcePath('/DbAmbigua/' + db)
        FinalResult = self.resourcePath('/FinalResult/' + db)
        if total:
            try:
                shutil.rmtree(Database)
            except:
                pass
        try:
            shutil.rmtree(BlastDb)
        except:
            pass
        try:
            shutil.rmtree(BlastResult)
        except:
            pass
        try:
            shutil.rmtree(DbAmbigua)
        except:
            pass
        try:
            shutil.rmtree(FinalResult)
        except:
            print("Error in a deletion")

    def deleteSeq(self, db, seqPath):
        try:
            os.remove(seqPath)
            self.restartDb(db)
        except:
            print("La carpeta no existe")
        self.configureDb(db)


    def addSeq(self, path,db, name, seq):
        file = open(path+"/"+name+".fa", "w")
        file.write(">"+name + os.linesep)
        file.write(seq)
        file.close()
        self.restartDb(db)


    def setDb(self, dbName):
        self.db = dbName
        self.simpleDbCreator = SC.SimpleDbCreator("Databases/"+self.db, "Blastdb",
                                                  self.db, "secuencias", "fasta")
        self.globalBlast = GC.GlobalBlast("Blastdb", "secuencias", "salida", "fasta", "BlastResult", self.db)
        self.ambiguousDbCreator = AC.AmbiguousDbCreator("BlastResult", "Nuevadb" , "salida", "fasta", "DbAmbigua", self.db)
        self.simpleBlast = S.SimpleBlast("DbAmbigua", "salida", "salida", "fasta", "FinalResult", self.db, True)
        self.resultsAnalizer = RA.ResultsAnalizer("FinalResult",self.db, True)

    def setAddSeqValues(self, path, content,name):
        self.path = path
        self.newSeqContent = content
        self.newSeqName = name

    def setSequenceToDelete(self, seqPath):
        self.sequenceToDelete = seqPath

    def restartDb(self,db):
        self.deleteDb(db, False)
        self.configureDb(db)
