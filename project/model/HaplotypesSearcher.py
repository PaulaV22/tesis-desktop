import SimpleDbCreator as SC
import SimpleBlast as S
from project.model import AmbiguousDbCreator as AC, GlobalBlast as GC
import ResultsAnalizer as RA
import os
import shutil
import json
from Signal import HaploSignal
import sys
import MailSender as Mail
from PySide.QtGui import *
from PySide.QtCore import *

class HaplotypesSearcher(QRunnable):

    def __init__(self, dbName=None):
        QRunnable.__init__(self)
        if (dbName):
            self.db=dbName
        else:
            self.db = "BoLa"
        self.projectPath = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

        self.categories  = {"ALTA": 1, "MEDIA": 1, "BAJA": 1}
        self.categoriesPath = self.projectPath+"/Categories"
        self.setDb(self.db)
        if not os.path.exists(self.categoriesPath):
            os.makedirs(self.categoriesPath)
        self.signals = HaploSignal()
        self.option="compare"
        self.newDb=None
        self.mailsvc = Mail.MailSender()
        #self.dbAdmin = DB.DbAdmin()

    def resourcePath(self,relative_path):
        """ Get absolute path to resource, works for dev and for PyInstaller """
        base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        output = base_path + relative_path
        return output

    def setOption(self,option):
        self.option=option

    def setNewDb(self,newDb):
        self.newDb=newDb

    def getResults(self,queryName,queryPath, database, numSeqs):
        simpleBlast = S.SimpleBlast("DbAmbigua", "salida", "salida", "fasta", "FinalResult", database, True)
        simpleBlast.align(queryPath, queryName)
        resultsAnalizer = RA.ResultsAnalizer("FinalResult",database, self.categories,self.categoriesPath)

        results = resultsAnalizer.getSimilarSequences(queryName,numSeqs)
        return results

    def run(self):
    #def searchHaplotypes(self):
        if self.option=="compare":
            #queryPath = self.projectPath + "/" +self.queryName
            queryPath = self.resourcePath("/"+self.queryName)
            # alinear y obtener resultados de la query deseada
            self.simpleBlast.align(queryPath, self.queryName)
            results = self.resultsAnalizer.getSimilarSequences(self.queryName,self.numSeqs)
            self.signals.result.emit(results)
            return results
            #return results
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

    def searchHaplotypes(self, queryName, numSeq):
        queryPath = self.projectPath + "/" + queryName

        # alinear y obtener resultados de la query deseada
        self.simpleBlast.align(queryPath, queryName)
        self.resultsAnalizer.getSimilarSequences(queryName, numSeq)
        print ("fin")

    def setAddSeqValues(self, path, content,name):
        self.path = path
        self.newSeqContent = content
        self.newSeqName = name

    def setSequenceToDelete(self, seqPath):
        self.sequenceToDelete = seqPath

    def setQueryData(self,queryName, numSeqs):
        self.queryName = queryName
        self.numSeqs = numSeqs

    def probarGlobalComparator(self):
        self.globalBlast.align("BoLa")

    def probarAmbiguousDbCreator(self):
        self.ambiguousDbCreator.makeDb()
        query = self.resourcePath('/BoLa/prueba.fa')
        self.simpleBlast.align(query, "prueba")

    def probarSimpleDbCreator(self):
        self.simpleDbCreator.makeDb()


    def deleteseqAdmin (self, db, sequence):
        self.dbAdmin.deleteSequence(self.projectPath,db,sequence)

    def configureDb(self, db):
        ####crear la bd con los archivos originales de BoLa####
        ready = False
        self.simpleDbCreator.makeDb()
        ####alinear todas las secuencias de BoLa entre si generando un archivo de salida por cada alineacion (n x n)####
        while not ready:
            try:
                self.globalBlast.align("/Databases/"+db)
                ready = True
            except:
                self.simpleDbCreator.makeDb()
                ready = False
        ####armar la base de datos con las posibles combinaciones (Nuevadb)####
        self.ambiguousDbCreator.makeDb()
        categories = {}
        #categoriesFile = self.projectPath + "/Categories/" + db + ".json"
        categoriesFile = self.resourcePath("/Categories/" + db + ".json")

        with open(categoriesFile, mode='w+') as f:
            json.dump(categories, f)

    def restartDb(self,db):
        self.deleteDb(db, False)
        self.configureDb(db)

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
            pass
        try:
            os.remove(self.resourcePath("/Categories/"+db+".json"))
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

    def setCategoryToFilesInDb(self, db, folder, category):
        folderPath = self.resourcePath("/Databases/"+db+"/"+folder)
        categoriesFile = self.resourcePath("/Categories/"+db+".json")
        with open(categoriesFile, mode='r') as json_file:
            data = json.load(json_file)
        for bases, dirs, files in os.walk(folderPath):
            for file in files:
                #newJson = '{"'+os.path.basename(file.name)+'" : "'+self.categories[category]+'"}'
                data[file[:-3]] = category
        with open(categoriesFile, 'w') as f:  # writing JSON object
            json.dump(data, f)

    def getDatabases(self):
        output = []
        dbs= self.resourcePath("\Databases")
        print(dbs)
        for dir in os.listdir(dbs):
            output.append(dir)
        output.sort()
        return output

    def getDb(self):
        return self.resourcePath('\Databases\\'+self.db)

    def getDbName(self):
        return self.db

    def setDb(self, dbName):
        self.db = dbName
        self.simpleDbCreator = SC.SimpleDbCreator("Databases/"+self.db, "Blastdb",
                                                  self.db, "secuencias", "fasta")
        self.globalBlast = GC.GlobalBlast("Blastdb", "secuencias", "salida", "fasta", "BlastResult", self.db)
        self.ambiguousDbCreator = AC.AmbiguousDbCreator("BlastResult", "Nuevadb" , "salida", "fasta", "DbAmbigua", self.db)
        self.simpleBlast = S.SimpleBlast("DbAmbigua", "salida", "salida", "fasta", "FinalResult", self.db, True)
        self.resultsAnalizer = RA.ResultsAnalizer("FinalResult",self.db, self.categories,self.categoriesPath)

    def getProjectPath(self):
        return self.projectPath

searcher = HaplotypesSearcher()
searcher.getDatabases()
#searcher.configureDb("BoLa")
#searcher.setCategoryToFilesInDb('BoLa', 'Mas_frecuentes', "ALTA")
#searcher.setCategoryToFilesInDb('BoLa', 'Menos_1%', "MEDIA")
#searcher.setCategoryToFilesInDb('BoLa', 'Menos_2%', "MEDIA")
#searcher.setCategoryToFilesInDb('BoLa', 'No_encontrados', "BAJA")

#searcher.searchHaplotypes()

#searcher.congifureDb()
#searcher.deleteSequence("BoLa", "DERB3_4501.fa")
#
#searcher.deleteSeq("BoLa","C:\Users\Paula\PycharmProjects\Haplotypes\BoLa\Mas_frecuentes\DRB3_2705.fa")
#searcher.addSeq("C:\Users\Paula\PycharmProjects\Haplotypes\BoLa\Mas_frecuentes","BoLa", "DRB3_2705",
                          #"GGAGTATTATAAGAGAGAGTGTCATTTCTTCAACGGGACCGAGCGGGTGCGGTTCCTGGACAGATGCTACACTAATGGAGAAGAGACCGTGCGCTTCGACAGCGACTGGGGCGAGTTCCGGGCGGTGACCGAGCTAGGGCGGCCGGACGCCGAGTACTGGAACAGCCAGAAGGACTTCCTGGAGGAGAGGCGGGCCGCGGTGGACAGGGTGTGCAGACACAACTACGGGGTCGTGGAGAGTTTCACTGTG")
#searcher.probarGlobalComparator()
#searcher.probarSimpleDbCreator()
#searcher.probarAmbiguousDbCreator()

