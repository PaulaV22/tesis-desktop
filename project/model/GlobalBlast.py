import SimpleBlast as SB
import os
# ESTA CLASE HACE UNA COMPARACION DE TODAS LAS SECUENCIAS DE LA BASE DE DATOS CONTRA TODAS GENERANDO UN ARCHIVO
# DE SALIDA CON NCBIBLASTCOMMANDLINE QUE ALINEA LAS SECUENCIAS PARA LUEGO GENERAR UNA NUEVA SECUENCIA CON LA COMBINACION
# DE ELLAS PARA GENERAR UNA NUEVA BASE DE DATOS DE POSIBLES SECUENCIAS POLIMORFICAS.
# RECIBE:
#       1: DBPATH ES LA CARPETA DONDE ESTAN TODOS LOS ARCHIVOS QUE DEBERA COMPARAR ENTRE SI (BLASTBD)
#       2: NEWDB ES LA CARPETA DONDE SE VAN A GUARDAR LOS ARCHIVOS CON LAS NUEVAS SECUENCIAS CREADA
#       3: OUTPUTFILE ES EL NOMBRE DEL ARCHIVO DE SALIDA QUE GENERARA CON EL RESUMEN DE TODAS LAS SECUENCIAS (SECUENCIAS)
#       4: OUTPUTFORMAT ES EL FORMATO DEL ARCHIVO DE SALIDA QUE TIENE TODAS LAS NUEVAS SECUENCIAS CREADAS (FASTA)

class GlobalBlast(SB.SimpleBlast):
    def __init__(self, dbPath, newDb, outputFile, outputFormat, outputPath, dbName):
        # recibe (BlastDb, secuencias, salida, fasta, BoLa)
        SB.SimpleBlast.__init__(self, dbPath, newDb, outputFile, outputFormat,outputPath, dbName)
        self.projectPath = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


    def align(self, queryDB):
        #queryDB es BoLa o la carpeta que contenga todas las secuencias documentadas originales
        # por cada secuencia de la base de datos comparar con la base de datos BlastDB
        queryDB = self.resourcePath("/"+queryDB)
        for bases, dirs, files in os.walk(queryDB):
            for file in files:
                queryName = file[:-3]
                sequence = bases + '/' + file
                SB.SimpleBlast.align(self,sequence,queryName)
