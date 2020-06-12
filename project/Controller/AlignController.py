from PySide.QtGui import *
from PySide.QtCore import *
import numpy as np
from Controller import Controller
import json
from project.model import Aligner as Aligner
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle
import subprocess
import os

class AlignController(Controller):

    def __init__(self,window):
        Controller.__init__(self,window)
        self.threadPool = QThreadPool()

    def configureView(self):
        self.window.buttonAlign.clicked.connect(self.align)
        self.window.progressBar_5.setMinimum(0)
        self.window.progressBar_5.setMaximum(0)
        self.window.progressBar_5.setVisible(False)
        self.window.buttonSeeAlignment.setVisible(False)
        self.window.buttonSeeAlignment.clicked.connect(self.downloadPDF)
        # self.HS.signals.result.__call__(self.showResults)
        self.HS.signals.updatedatabases.connect(self.dbReady)
        self.HS.signals.aligned.connect(self.showAlignment)

    def setSequences(self, seq1Name, seq1Content, seq2Name, seq2Content):
        self.db = "align"
        self.seq1Name = seq1Name
        self.seq1Content = seq1Content
        self.seq2Name = seq2Name
        self.seq2Content = seq2Content


    def align(self):
        self.window.buttonAlign.hide()
        self.window.progressBar_5.show()
        self.window.buttonAlign.repaint()
        self.window.progressBar_5.repaint()
        inputSequence1 = self.window.sequenceInput1.toPlainText()
        inputSequence2 = self.window.sequenceInput2.toPlainText()

        seq1Name = inputSequence1.partition('\n')[0]
        seq1Name = ''.join(e for e in seq1Name if e.isalnum())
        seq1Name = seq1Name + ".fa"
        self.aligner = Aligner.Aligner()
        self.aligner.setProjectPath(self.getProjectPath())
        self.aligner.setSequences(seq1Name, inputSequence1, inputSequence2)
        self.aligner.signals.aligned.connect(self.showAlignment)

        tempFile = self.HS.getProjectPath()+"/tmp.fa"
        file = open(tempFile, "w+")
        file.write(inputSequence1)
        file.close()
        self.threadPool.start(self.aligner)
        #self.HS.searchHaplotypes("tmp.fa", self.numSeqs)


    def showAlignment(self, results):
        self.window.progressBar_5.setVisible(False)
        self.results = results
        self.showResults()
        self.window.buttonSeeAlignment.setVisible(True)





    def dbReady(self):
        self.setDatabases()

    def downloadPDF(self):
        print("Download")


    def getContent(self, sequence):
        arr = sequence.split("\n")
        arr.pop(0)
        content = ''.join(arr)
        return content

    def min(self, val1, val2):
        if (val1 >= val2): return val2
        return val1

    def max(self, val1, val2):
        if (val1 >= val2): return val1
        return val2


    def showResults(self):
        self.rows1 = []
        self.rows2 = []
        if (len(self.results) > 0):
            results = json.loads(self.results)
            alignment = results[len(results)-1]

            #FIRST TABLE
            score =  alignment['score']
            evalue = alignment['evalue']
            similarity = alignment['similarity']

            r1 =["SCORE", "EVALUE", "SIMILARITY"]
            r2 = [score, evalue, similarity]
            self.rows1.append(r1)
            self.rows1.append(r2)


            array = alignment['alignment'].split('>')
            sequence1 = array[1]
            sequence2 = array[2]
            sequence1Name = sequence1.split('\n')[0]
            sequence1Name = sequence1Name.split(' ')[0]

            sequence2Name = sequence2.split('\n')[0]
            sequence2Name = sequence2Name.split(' ')[0]

            sequence1Content = self.getContent(sequence1)
            sequence2Content = self.getContent(sequence2)

            seq1start = alignment['queryStart']
            seq2start = alignment['hitStart']
            length = 50
            iteration = 0

            middleSeq = ""
            max = self.max(len(sequence1Content), len(sequence2Content))
            for j in range (0,max-1):

                if (sequence1Content[j] and sequence2Content[j]):
                    if (sequence1Content[j] == sequence2Content[j]):
                        middleSeq = middleSeq + "|"
                    else:
                        middleSeq = middleSeq + "."
                else:
                    middleSeq = middleSeq + "."


            #DIVIDIR EN ARREGLO DE STRINGS DE 50 ELEMENTOS (VER SI FUNCIONA)
            chunk_size = 50
            chunkedSeq1 =[sequence1Content[i:i + chunk_size] for i in range(0, len(sequence1Content), chunk_size)]
            chunkedSeq2 =[sequence2Content[i:i + chunk_size] for i in range(0, len(sequence2Content), chunk_size)]
            chunkedMiddle =[middleSeq[i:i + chunk_size] for i in range(0, len(middleSeq), chunk_size)]


            for i in range((len(chunkedMiddle))):

                s1contentValue = " "
                s1StartValue = " "
                s1endValue = " "
                if (chunkedSeq1[i]):
                    s1contentValue= chunkedSeq1[i]
                    s1StartValue = seq1start + len(chunkedSeq1[i]) * i+1
                    s1endValue = s1StartValue-1 + len(chunkedSeq1[i])



                row1 = []
                row1.append(sequence1Name)
                row1.append(s1StartValue)
                row1.append(s1contentValue)
                row1.append(s1endValue)


                row2 = []
                row2.append(" ")
                row2.append(" ")
                row2.append(chunkedMiddle[i])
                row2.append(" ")

                s2StartValue = " "
                s2contentValue = " "
                s2endValue = " "
                if (chunkedSeq2[i]):
                    s2contentValue = chunkedSeq2[i]
                    s2StartValue = seq2start + len(chunkedSeq2[i]) * i + 1
                    s2endValue = s2StartValue - 1 + len(chunkedSeq2[i])


                row3 = []
                row3.append(sequence2Name)
                row3.append(s2StartValue)
                row3.append(s2contentValue)
                row3.append(s2endValue)

                self.rows2.append(row1)
                self.rows2.append(row2)
                self.rows2.append(row3)

        print(self.rows1)
        print(self.rows2)

    def downloadPDF(self):
        if os.path.exists("output.pdf"):
            os.remove("output.pdf")
        doc = SimpleDocTemplate("output.pdf", pagesize=letter)
        # container for the 'Flowable' objects
        elements = []
        data= [['00', '01', '02', '03', '04'],
        ['10', '11', '12', '13', '14'],
        ['20', '21', '22', '23', '24'],
        ['30', '31', '32', '33', '34']]
        t1 = Table(self.rows1)
        t1.setStyle(TableStyle([(('FONTNAME', (0, 0), (-1, -1), 'Courier'))]))
        elements.append(t1)

        t2=Table(self.rows2)
        t2.setStyle(TableStyle([(('FONTNAME', (0,0), (-1,-1), 'Courier'))]))
        elements.append(t2)
        # write the document to disk
        doc.build(elements)
        subprocess.Popen(["output.pdf"], shell=True)