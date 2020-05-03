from __future__ import division
import os
from Bio import SearchIO
import sys
import json

class ResultsAnalizer():

    def __init__(self, resultsPath, dbName, categories, categoriesPath):
        # recibe (FinalResult, BoLa o Prueba2)
      #  self.projectPath = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.resultsPath = self.resourcePath("/"+resultsPath)
        self.dbName = dbName
        self.resultFiles = self.resourcePath("/"+resultsPath+"/"+dbName)
        self.categories = categories
        self.categoriesPath = categoriesPath

    def resourcePath(self,relative_path):
        """ Get absolute path to resource, works for dev and for PyInstaller """
        base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        output = base_path + relative_path
        return output

    def mini(self, weight1, weight2):
        if (weight1>weight2):
            return weight2
        return weight1

    def getWeight(self, sequenceId, sequencesCategories, file):
        seq1 = sequenceId.split('-')[0]
        seq2 = sequenceId.split('-')[1]
        seq1Id = seq1.replace("*", "_")
        seq2Id = seq2.replace("*", "_")
        if ("DRB3*902/DRB3*1103" in sequenceId):
            weight1 = sequencesCategories.get(seq1Id)
            weight2 = sequencesCategories.get(seq2Id)
            if (weight1 is None) and (weight2 is None):
                return 1
            if not (weight1 is None) and (weight2 is None):
                return self.categories.get(weight1)
            if not (weight2 is None) and (weight1 is None):
                return self.categories(weight2)
            if not (weight1 is None) and not (weight2 is None):
                return self.mini(self.categories.get(weight1), self.categories.get(weight2))
        else:
            weight1 = sequencesCategories.get(seq1Id)
            weight2 = sequencesCategories.get(seq2Id)
            if (weight1 is None) and (weight2 is None):
                return 1
            if not (weight1 is None) and (weight2 is None):
                return self.categories.get(weight1)
            if not (weight2 is None) and (weight1 is None):
                return self.categories(weight2)
            if not (weight1 is None) and not (weight2 is None):
                return self.mini(self.categories.get(weight1), self.categories.get(weight2))

    def getStandarName(self,id):
        id = id.replace("_", "*")
        id = id.replace("/", "-")
        return id

    def getComplementary(self, id):
        seq1 = id.split('-')[0]
        seq2 = id.split('-')[1]
        output = seq2+"-"+seq1
        return output

    def getSimilarSequences(self, name, number=None):
        self.resultFiles = str(self.resultFiles) + "/"+str(name)
        i = 0
        sequences = dict()
        #categoriesFile = self.categoriesPath+"/"+self.dbName +".json"
        #with open(categoriesFile, mode='r') as json_file:
        #    sequencesCategories = json.load(json_file)
        for bases, dirs, files in os.walk(self.resultFiles):
            for file in files:
                outputName = bases + "/" + file
                result = SearchIO.read(outputName, "blast-xml")
                for hits in result:
                    i = i+1
                    hsp = hits[0]
                    id = hsp.hit.id
                    score = hsp.bitscore
                    evalue = hsp.evalue
                    positives = hsp.pos_num
                    align_length = hsp.aln_span
                    percent = float("{0:.3f}".format(positives/align_length))
                    id = self.getStandarName(id)
                    #weight = self.getWeight(id,sequencesCategories, file)
                    complementary = self.getComplementary(id)
                    if not (complementary in sequences.keys()):
                        #sequences[id] = [score * weight, evalue, positives, align_length, percent*weight]
                        sequences[id] = [score, evalue, positives, align_length, percent]

        n = 0
        salida = []
        for key, value in sorted(sequences.items(), key=lambda item: item[1][4], reverse=True):
            if (n<int(number)):
                value.insert(0, key)
                print(value)
                data = {'id' :value[0], 'score':value[1], 'evalue':value[2], 'similarity':value[5]}
                dataJson = json.dumps(data)
                salida.append(data)
            n = n +1
        #print(salida)
        salidaJson = json.dumps(salida)
        #print(number)
        return salidaJson

