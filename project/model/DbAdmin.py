import os
from Bio import SeqIO
from Bio import SearchIO
import shutil
from shutil import copyfile
import subprocess
import SimpleBlast as SB

class DbAdmin:

    def deleteSequence(self, projectPath, db, file):
        self.projectPath = projectPath
        self.db = db
        self.deleteFromFolder(projectPath,db,file)
        self.deleteFromDb(projectPath,db,file, "BlastDb")
        self.deleteFromAlignResult(projectPath,db,file,"BlastResult")
        self.deleteFromDb(projectPath, db, file, "DbAmbigua")

    def contains(self, string1, string2):
        s1=(string1.replace('_', ''))
        s1=(s1.replace('*',''))
        s2=(string2.replace('_', ''))
        s2=(s2.replace('*',''))
        return s2 in s1

    # Borrarlo de la carpeta contenedora
    def deleteFromFolder(self,projectPath, db, file):
        dbPath = projectPath+"/"+db
        for bases, dirs, files in os.walk(dbPath):
            for f in os.listdir(bases):
                filePath = bases + '/' + f
                if self.contains(f, file):
                    try:
                        os.remove(filePath)
                        return
                    except Exception as e:
                        print e

    # Borrarlo de la base de datos BlastDb o AmbiguousDb
    def deleteFromDb(self,projectPath, db,file, dbType ):
        dbPath = projectPath + "/" + dbType+"/" + db
        for bases, dirs, files in os.walk(dbPath):
            for f in os.listdir(bases):
                filePath = bases + '/' + f
                # si es un archivo fasta y no un directorio  --> ver, creo que hay una manera mas elegante de preguntar si es un archivo o directorio
                if self.contains(filePath, file[:-3]):
                    if os.path.isfile(filePath):
                        os.remove(filePath)
                    else:
                        shutil.rmtree(filePath)
                else:
                    if os.path.isfile(filePath) and self.contains(f,"fasta"):
                        with open(filePath) as originalFasta, open(bases+ "/output.fasta", 'w') as correctedFasta:
                            for seq_record in SeqIO.parse(originalFasta, "fasta"):
                                if not(self.contains(seq_record.id,file[:-3])):
                                    SeqIO.write(seq_record, correctedFasta, 'fasta')
                        os.remove(filePath)
                        try:
                            os.rename(bases+"/output.fasta", filePath)
                        except Exception as e:
                            print(e)
                        folder = os.path.dirname(filePath)
                        self.removeOtherFilesInFolder(folder)
                        filename= os.path.splitext(f)[0]
                        output = folder+"/"+filename
                        dbName = os.path.basename(folder)
                        while (self.testFails(folder, dbName, filename)):
                            self.generateBlastDb(filePath, output)


    def removeOtherFilesInFolder(self,folder):
        for file in os.listdir(folder):
            if not self.contains(file,".fasta"):
                os.remove(folder+"/"+file)

    def generateBlastDb(self, filePath, output):
        command = 'powershell.exe makeblastdb -in ' + filePath + ' -out ' + output + ' -parse_seqids -dbtype nucl'
        subprocess.Popen(command)
        print (subprocess.check_output(command))

    def testFails(self,db):
        print("A TESTEAR DB "+db)
        if (len([f for f in os.listdir(db)]) <7):
            return True
        return False

    def testFails(self,db, dbName,file):
        print("A TESTEAR DB "+db)
        if (len([f for f in os.listdir(db)]) <7):
            return True
        outputPath = self.projectPath+"/Test"+"/"+dbName+"/"+file
        sb = SB.SimpleBlast(dbName, "salida", "salida", "fasta", "Test", self.db)
        queryName = "queryTestBola.fa"
        queryPath = self.projectPath + "/" + queryName
        try:
            sb.align(queryPath, queryName)
        except Exception as e:
            print(e)
            return True
        return False

    # Borrarlo de la base de datos ambigua BlastResult
    def deleteFromAlignResult(self, projectPath, db, file, alignResult):
        blastResultPath = projectPath + "/" + alignResult + "/" + db
        for bases, dirs, files in os.walk(blastResultPath):
            for f in os.listdir(bases):
                filePath = bases + '/' + f
                # si el nombre del archivo a leer contiene el nombre del archivo a borrar, borrarlo
                if (self.contains(f, file[:-3])):
                    os.remove(filePath)
                # sino puede estar adentro la alineacion
                else:
                    with open(filePath) as originalXml:
                        result = SearchIO.read(originalXml, "blast-xml")
                        i = 0
                        for hits in result:
                            hsp = hits[0]
                            id = hsp.hit.id
                            if (id == file[:-3]):
                                result.hit_keys.pop(i)
                                result.hits.pop(i)
                                result.hsps.pop(i)
                                result.pop(i)
                            i = i + 1
                    os.remove(filePath)
                    SearchIO.write(result, filePath, 'blast-xml')

    def addSequence(self, source, projectPath, db, subPath=None, ):
        file = os.path.basename(source)
        self.addToFolder(source, projectPath,db,file, subPath)


    def addToFolder(self, source,projectPath, db,file,subpath):
        if  not subpath is None:
            fullPath = projectPath+"/"+db+"/"+subpath+"/"+file
        else:
            fullPath = projectPath+"/"+db+"/"+file
        copyfile(source, fullPath)

   # def addToDatabase(self, source,projectPath, db, file,subpath, dbType):


