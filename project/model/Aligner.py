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
import SimpleDbCreator as SC
import HaplotypesSearcher as HaplotypesSearcher
import time
from project.model import DbAdmin as DbAdmin

class Aligner(QRunnable):

    def __init__(self, dbName=None):
        QRunnable.__init__(self)
        self.HS = HaplotypesSearcher.HaplotypesSearcher()
        self.signals = HaploSignal()

    def setProjectPath(self, projectPath):
        self.projectPath = projectPath

    def setSequences(self, seq1Name, seq1Content, seq2Content):
        self.db = "align"
        self.seq1Name = seq1Name
        self.seq1Content = seq1Content
        self.seq2Content = seq2Content

    def run(self):
        dbPath = self.projectPath + "/Databases/" + self.db
        if os.path.exists(dbPath):
            shutil.rmtree(dbPath)
        os.makedirs(dbPath)
        seq1Path = dbPath + "/" + self.seq1Name
        while not os.path.exists(dbPath):
            time.sleep(1)

        file = open(seq1Path, 'w+')
        file.write(self.seq1Content)
        file.close()

        simpleDbCreator = SC.SimpleDbCreator("Databases/" + self.db, "Blastdb", self.db, "secuencias", "fasta")
        simpleDbCreator.makeDb()

        results = self.align()
        self.signals.aligned.emit(results)
        self.dbAdmin = DbAdmin.DbAdmin(self.db)
        self.dbAdmin.deleteDb(self.db)
        self.signals.alignedDeleted.emit()
        return results


    def align(self):
        if not self.HS:
            self.HS = HaplotypesSearcher.HaplotypesSearcher()
        self.HS.setDb(self.db)
        seqName = self.seq2Content.partition('\n')[0]
        seqName = ''.join(e for e in seqName if e.isalnum())

        #seqPath = "/Temporal/" + seqName + ".fa"
        #tempFile = self.projectPath + seqPath

        #print("TEMP FILE IN PATH " + tempFile)
        seqPath = "/tmp.fa"
        tempFile = self.projectPath + seqPath

        print("TEMP FILE IN PATH " + tempFile)
        file = open(tempFile, 'w+')
        file.write(self.seq2Content)
        file.close()

        results = self.HS.getResults(seqName, tempFile, self.db, 1, False)
        try:
            os.remove(tempFile)
        except Exception as e:
            print(e)
        return results

