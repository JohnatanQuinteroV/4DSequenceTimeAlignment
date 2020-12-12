#Step4.py "4D sequence time alignment" 
#Johnatan Quintero V. johnatan2425@hotmail.com

#Step4: Toma dos archivos de texto textA.txt y textB.txt, con los codigos de A y B
#       Utiliza una libreria eficiente para secuencias largas, Align de BioPython  
#       Se genera 1 archivo con el resultado de la alineacion, textC.txt

from Bio import Align  #libreria hecha en C para alinear secuencias largas
import numpy as np   
import pandas as pd     

#lee strings de codigos de textA y textB
#retorna la alineacion en un texto llamado textC.txt
def alignment():
    file = open('textA.txt', 'r')  #leyendo string de sec A
    sa, ssa = '', ''
    for line in file.readlines():
        sa = line
        ssa = ssa+sa  #string de sec A
    file.close()
    
    file = open('textB.txt', 'r')  #leyendo string de sec B
    sb, ssb = '', ''
    for line in file.readlines():
        sb = line
        ssb = ssb+sb  #string de sec B
    file.close()
    
    aligner = Align.PairwiseAligner()  
    alignments = aligner.align(ssa, ssb)  #alineando strings
    aline = alignments[0]  #solo guardamos 1 alineacion con el mejor score, ahorro de memoria
    file = open('textC.txt', 'w')#guardando string de la alineacion en textC.txt
    file.write(str(aline))
    file.close()

def main():
    alignment()
main =main()
