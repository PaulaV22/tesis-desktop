import os
from project.model import HaplotypesSearcher as HaplotypeSearcher
import sys

class Controller():
    def __init__(self, window):
        self.window = window
        self.projectPath = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.HS = HaplotypeSearcher.HaplotypesSearcher()
        self.dbList = self.getDatabases()
        self.dbName =""

    def resourcePath(self,relative_path):
        """ Get absolute path to resource, works for dev and for PyInstaller """
        base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        output = base_path + relative_path
        return output

    def getProjectPath(self):
        return self.projectPath

    def getDatabases(self):
        output = []
        dbs= self.resourcePath("/Databases")
        for dir in os.listdir(dbs):
            output.append(dir)
        output.sort()
        return output


    def getDb(self):
        self.dbList = self.getDatabases()
        selectedIndex = self.window.selectDatabase.currentIndex()
        dbName = self.dbList[selectedIndex]
        return dbName


    def setDb(self, dbName):
        self.dbName= dbName
        self.HS.setDb(dbName)


    def setDatabases(self):
        self.dbList = self.HS.getDatabases()
        self.window.selectDeleteDatabase.clear()
        self.window.selectDatabase.clear()
        for db in self.dbList:
            self.window.selectDatabase.addItem(db)
            self.window.selectDeleteDatabase.addItem(db)