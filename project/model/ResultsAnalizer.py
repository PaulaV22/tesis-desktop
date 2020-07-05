from __future__ import division
import os
from Bio import SearchIO
import sys
import json
import shutil
class ResultsAnalizer():

    def __init__(self, resultsPath, dbName, ambiguo = True):
        # recibe (FinalResult, BoLa o Prueba2)
      #  self.projectPath = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.resultsPath = self.resourcePath("/"+resultsPath)
        self.dbName = dbName
        self.resultFiles = self.resourcePath("/"+resultsPath+"/"+dbName)
        self.ambiguo = ambiguo

    def resourcePath(self,relative_path):
        """ Get absolute path to resource, works for dev and for PyInstaller """
        base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        output = base_path + relative_path
        return output

    def mini(self, weight1, weight2):
        if (weight1>weight2):
            return weight2
        return weight1


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
                    if not hasattr(hsp, 'pos_num'):
                        percent = 1
                    else:
                        positives = hsp.pos_num
                        align_length = hsp.aln_span
                        percent = float("{0:.3f}".format(positives/align_length))
                    complementary = ""
                    if self.ambiguo:
                        id = self.getStandarName(id)
                        complementary = self.getComplementary(id)
                    alignment = hsp.fragment.aln
                    if not (complementary in sequences.keys()):
                        # sequences[id] = [score * weight, evalue, positives, align_length, percent*weight]
                        sequences[id] = [score, evalue, align_length, percent, alignment.format("fasta"),
                                         hsp.query_start, hsp.hit_start]
        n = 0
        salida = []
        #for key, value in sorted(sequences.items(), key=lambda item: item[1][4], reverse=True):
        for key, value in sorted(sequences.items(), key=lambda item: item[1][3], reverse=True):
            if (n<int(number)):
                value.insert(0, key)
                print(value)
                data = {'id': value[0], 'score': value[1], 'evalue': value[2], 'similarity': value[4],
                        'alignment': value[5], 'queryStart': value[6], 'hitStart': value[7]}
                salida.append(data)
            n = n +1
        #print(salida)
        salidaJson = json.dumps(salida)
        #print(number)
        shutil.rmtree(self.resultsPath)
        return salidaJson
