from PySide.QtGui import *
from PySide.QtCore import *
import numpy as np
from Controller import Controller
import json
from project.model import HaplotypesSearcher as HaplotypeSearcher

class CustomTableModel(QAbstractTableModel):
    def __init__(self, data=None):
        QAbstractTableModel.__init__(self)
        self.input_seqs = list()
        self.input_score = list()
        self.input_evalue = list()
        self.input_percent = list()
        self.load_data(data)


    def load_data(self, data):
        if len(data)>1:
            for dic in data:
                self.input_seqs.append(dic['id'])
                self.input_score.append(dic['score'])
                self.input_evalue.append(dic['evalue'])
                self.input_percent.append(dic['similarity'])
                self.column_count = 4
                self.row_count = len(self.input_seqs)
        else:
            self.column_count = 4
            self.row_count = 0

    def rowCount(self, parent=QModelIndex()):
        return self.row_count

    def columnCount(self, parent=QModelIndex()):
        return self.column_count

    def headerData(self, section, orientation, role):
        if role != Qt.DisplayRole:
            return None
        if orientation == Qt.Horizontal:
            return ("Combinacion", "Score", "E-Value","Similitud")[section]
        else:
            return "{}".format(section)

    def data(self, index, role = Qt.DisplayRole):
        column = index.column()
        row = index.row()
        if role == Qt.DisplayRole:
            if column == 0:
                return self.input_seqs[row]
            elif column == 1:
                return self.input_score[row]
            elif column == 2:
                return self.input_evalue[row]
            elif column == 3:
                return self.input_percent[row]
        elif role == Qt.BackgroundRole:
            return QColor(Qt.white)
        elif role == Qt.TextAlignmentRole:
            return Qt.AlignRight
        return None

class CompareController(Controller):

    def __init__(self,window):
        Controller.__init__(self,window)
        self.threadPool = QThreadPool()

    def configureView(self):
        self.drawTable([[]])
        self.window.buttonCompare.clicked.connect(self.compare)
        self.window.resultLenghtSpin.setValue(10)
        self.window.progressBar.setMinimum(0)
        self.window.progressBar.setMaximum(0)
        self.window.progressBar.setVisible(False)
        # self.HS.signals.result.__call__(self.showResults)
        self.HS.signals.updatedatabases.connect(self.dbReady)
        self.HS.signals.result.connect(self.showResults)


    def compare(self):
        self.drawTable([[]])
        self.window.buttonCompare.hide()
        self.window.progressBar.show()
        self.window.buttonCompare.repaint()
        self.window.progressBar.repaint()
        inputSequence = self.window.sequenceInput.toPlainText()
        self.HS = HaplotypesSearcher.HaplotypesSearcher(self.dbName)
        self.HS.signals.result.connect(self.showResults)
        tempFile = self.HS.getProjectPath()+"/tmp.fa"
        file = open(tempFile, "w+")
        file.write(inputSequence)
        file.close()
        self.numSeqs = self.window.resultLenghtSpin.value()
        self.HS.setQueryData("tmp.fa", self.numSeqs)
        self.HS.setDb(self.dbName)
        self.HS.setOption("compare")
        self.threadPool.start(self.HS)
        #self.HS.searchHaplotypes("tmp.fa", self.numSeqs)

    def drawTable(self, data):
        self.tablemodel = CustomTableModel(data)
        self.window.tableResults.setModel(self.tablemodel)
        self.horizontal_header = self.window.tableResults.horizontalHeader()
        self.horizontal_header.resizeSections(QHeaderView.Stretch)
        self.horizontal_header.setStretchLastSection(True)


    def showResults(self, results):
        print(results)
        #dt = np.dtype('string', 'float', 'int')

        #resultsArray = np.array(results, dtype=dt)
        resultsArray = json.loads(results)
        print(resultsArray)
        model = CustomTableModel(resultsArray)
        # Creating a QTableView
        self.window.tableResults.setModel(model)
        self.window.progressBar.hide()
        self.window.buttonCompare.show()
        # self.drawTable(results)

    def compare(self):
        self.drawTable([[]])
        self.window.buttonCompare.hide()
        self.window.progressBar.show()
        self.window.buttonCompare.repaint()
        self.window.progressBar.repaint()
        self.dbName = self.getDb()
        inputSequence = self.window.sequenceInput.toPlainText()
        self.HS = HaplotypeSearcher.HaplotypesSearcher(self.dbName)
        self.HS.setDb(self.dbName)
        self.HS.signals.result.connect(self.showResults)
        tempFile = self.HS.getProjectPath()+"/tmp.fa"
        file = open(tempFile, "w+")
        file.write(inputSequence)
        file.close()
        self.numSeqs = self.window.resultLenghtSpin.value()
        self.HS.setQueryData("tmp.fa", self.numSeqs)
        self.HS.setOption("compare")
        self.threadPool.start(self.HS)
        #self.HS.searchHaplotypes('tmp', self.numSeqs)
        #results = [['DRB3*0201-DRB3*2704', 235.645, 1.57613e-65, 140, 147, 0.952], ['DRB3*1201-DRB3*2704', 239.338, 1.22034e-66, 143, 151, 0.947], ['DRB3*1201-DRB3*0201', 239.338, 1.22034e-66, 139, 147, 0.946], ['DRB3*0201-DRB3*0201', 239.338, 1.21842e-66, 139, 147, 0.946], ['DRB3*701-DRB3*0201', 228.258, 2.63927e-63, 139, 147, 0.946], ['DRB3*0201-DRB3*2708', 235.645, 1.57613e-65, 139, 147, 0.946], ['DRB3*0201-DRB3*2701', 235.645, 1.57613e-65, 139, 147, 0.946], ['DRB3*0201-DRB3*2710', 235.645, 1.57613e-65, 139, 147, 0.946], ['DRB3*1201-DRB3*1201', 241.185, 3.393e-67, 142, 151, 0.94], ['DRB3*701-DRB3*1201', 231.952, 2.04028e-64, 142, 151, 0.94]]

    def drawTable(self, data):
        #header = ["Combinacion", "Score", "E-value", "Similitud"]
        self.tablemodel = CustomTableModel(data)

        # Creating a QTableView
        self.window.tableResults.setModel(self.tablemodel)
        self.horizontal_header = self.window.tableResults.horizontalHeader()
        self.horizontal_header.resizeSections(QHeaderView.Stretch)
        self.horizontal_header.setStretchLastSection(True)

    def dbReady(self):
        self.setDatabases()
