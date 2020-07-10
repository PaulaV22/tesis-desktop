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

    def getResults(self,queryName,queryPath, database, numSeqs, ambiguo = True):
        if ambiguo:
            db = "DbAmbigua"
        else:
            db = "Blastdb"
        simpleBlast = S.SimpleBlast(db, "salida", "salida", "fasta", "FinalResult", database, True)
        simpleBlast.align(queryPath, queryName)
        resultsAnalizer = RA.ResultsAnalizer("FinalResult", database, ambiguo)
        results = resultsAnalizer.getSimilarSequences(queryName, numSeqs)
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


    def searchHaplotypes(self, queryName, numSeq):
        queryPath = self.projectPath + "/" + queryName

        # alinear y obtener resultados de la query deseada
        self.simpleBlast.align(queryPath, queryName)
        self.resultsAnalizer.getSimilarSequences(queryName, numSeq)
        print ("fin")

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
        self.resultsAnalizer = RA.ResultsAnalizer("FinalResult",self.db, True)

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

