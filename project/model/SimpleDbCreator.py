from Bio import SeqIO
import os
import subprocess
import sys

from project.model import DbCreator as DBC

# ESTA CLASE CREA UNA BASE DE DATOS BLAST A PARTIR DE ARCHIVOS DE SECUENCIAS.
# PARA ESO RECIBE:
#       1: EL PATH DONDE ESTAN LOS ARCHIVOS FASTA DE CADA SECUENCIA
#       2: EL NOMBRE DE LA NUEVA BASE DE DATOS BLAST A CREAR (Blastdb)
#       3: EL NOMBRE DE LA BASE DE DATOS SELECCIONADA (BOLA)
#       3: OUTPUTFILE ES EL NOMBRE DE UN ARCHIVO INTERMEDIO QUE REUNE TODAS LAS SECUENCIAS (CAMBIAR MAS REPRESENTATIVO)
#       4: OUTPUT FORMAT ES FASTA


class SimpleDbCreator(DBC.DbCreator):
    def __init__(self, filesPath, newDb, dbName, outputFile, outputFormat):
        DBC.DbCreator.__init__(self, filesPath, newDb, dbName, outputFile, outputFormat)



        #self.projectPath = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        # directorio donde se encuentran las secuencias individuales
        self.filesPath = self.resourcePath("/"+ filesPath)
        # directorio donde se creara la base de datos
        self.newDb = newDb
        #Nommbre de la base de datos seleccionada
        self.dbName = dbName
        # archivo intermedio usado para crear la base de datos
        self.outputFile = outputFile
        self.outputFormat = outputFormat

    def resourcePath(self,relative_path):
        """ Get absolute path to resource, works for dev and for PyInstaller """
        base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        output = base_path + relative_path
        return output

    def createFolder(self, newFolder):
        if not os.path.exists(newFolder):
            print ("creo el directorio "+ newFolder)
            os.makedirs(newFolder)

    def getSequences(self, filePath, sequences):
        for seq_record in SeqIO.parse(filePath, self.outputFormat):
            seq = seq_record.seq
            sequences.append(seq_record)  # o seq?
        return sequences


    def saveSequencesInFile(self, path, sequences, outputFile= None):
        if outputFile is None:
            with open(path + "/" + self.outputFile + "." + self.outputFormat, "w") as output_handle:
                SeqIO.write(sequences, output_handle, self.outputFormat)
        else:
            with open(path + "/" + outputFile + "." + self.outputFormat, "w") as output_handle:
                SeqIO.write(sequences, output_handle, self.outputFormat)

    def makeBlastDb(self, directory):
        # formo "secuencias.fasta"
        outputName = self.outputFile + "." + self.outputFormat
        if outputName in os.listdir(directory):
            dbpath1 = directory + '/' + outputName
            dbpath2 = directory + '/' + self.outputFile
            dbpath1.replace(" ","_")
            dbpath2.replace(" ","_")

            #print 'secuencias.fasta esta en ' + dbpath1
            # ver este comando que es el que tiene problemas
            #print (dbpath1)
            #print (dbpath2)
            command = 'powershell.exe makeblastdb -in ' + dbpath1 + ' -out ' + dbpath2 + ' -parse_seqids -dbtype nucl'
            subprocess.Popen(command)
            print (subprocess.check_output(command))



    def setFilesPath (self, filesPath):
        self.filesPath = filesPath

    def setOutputFile(self, file):
        self.outputFile = file

    def setOutputFormat(self, format):
        self.outputFormat= format

    def makeDb(self):
        dbPath = self.resourcePath('/'+self.newDb)
        dbPath= dbPath.replace(" ","_")
        self.createFolder(dbPath)
        # recorro los archivos y guardo las secuencias en un arreglo. Para cada directorio de las secuencias documentadas
        # creo el archivo secuencias.fasta para poder crear la base de datos en cada subdirectorio
        #print ("busca en los archivos " +self.filesPath)
        ##print (self.newDb)
        for bases, dirs, files in os.walk(self.filesPath):
            subdirectory = bases
            sequences = []
            bases = subdirectory.split('/')[2]
            newSubFolder = self.resourcePath("/"+self.newDb + "/" + bases)
            newSubFolder = newSubFolder.replace(" ", "_")

            # creo una subcarpeta para el subdirectorio correspondiente
            self.createFolder(newSubFolder)

            #print 'subfolder es ' + newSubFolder
            for file in os.listdir(subdirectory):
                filePath = subdirectory + '/' + file
                # si es un archivo fasta y no un directorio  --> ver, creo que hay una manera mas elegante de preguntar si es un archivo o directorio
                if os.path.isfile(filePath):
                    # obtengo las secuencias de cada archivo y la guardo en un arreglo de secuencias
                    sequences = self.getSequences(filePath, sequences)
                    # creo un archivo con las secuencias de ese subdirectorio
                    self.saveSequencesInFile(newSubFolder, sequences)
                    # con el archivo anterior puedo armar la base de datos en el subdirectorio
            self.makeBlastDb(newSubFolder)
