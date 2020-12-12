#Step2.py "4D sequence time alignment" 
#Johnatan Quintero V. johnatan2425@hotmail.com

#Step2: Toma archivos A.bvh, B.bvh y A_rot.csv, B_rot.csv. 
#       Usando matrices de rotacion se obtiene las posiciones locales.
#       Se generan 2 archivos, localesA.csv y localesB.csv 

import numpy as np    
import pandas as pd  
import csv         
import scipy as sc    
from scipy.spatial.transform import Rotation as R  #matrices de rotacion
              
#recibe archivo de secuencia X.bvh
#retorna arreglo con posiciones de offsets              
def offsets(archivo):      
    with open(archivo) as arch:   #se abre el archivo de la secuencia
        line = arch.readlines()   
        cont1, cont2 = 0, 0
        arrx = []
        for i in line:#recorriendo renglones 
            cont1 += 1     
            if 'ROOT' in i:#si esta la palabra ROOT
                a = i.strip()             #quita espacio en blanco inicial 
                a = a[4:len(a)]           #quita la palabra ROOT, deja el nombre
                b = line[cont1+1].strip() #baja dos renglones, quita el espacio
                b = b[6:len(b)]           #quita la palabra OFFSET, deja posiciones
                arrx.insert(cont2, a)     #guarda a, en lista arrx  
                cont2 += 1
                arrx.insert(cont2, b)     #guarda b, en lista arrx
                cont2 += 1
            elif 'JOINT' in i:#si esta la palabra JOINT
                c = i.strip()
                c = c[5:len(c)]           #quita la palabra JOINT..
                d = line[cont1+1].strip()
                d = d[6:len(d)]           #quita la palabra OFFSET..
                arrx.insert(cont2, c)     #guarda c, en lista arrx
                cont2 += 1
                arrx.insert(cont2, d)     #guarda d, en lista arrx
                cont2 += 1
            elif 'MOTION' in i:#si esta la palabra MOTION, sale del ciclo
                break
        arrayx = np.array(arrx)
        tmp=len(arrayx)/2
        posi = np.array_split(arrayx,tmp)   #modifica el arreglo en grupos de 3
        return posi      #retorna un arreglo con posiciones del offset
       
        
#recibe archivo de X_rot.csv de rotaciones. 
#retorna arreglo de rotaciones y lista con headers de marcadores. 
def rotaciones(archivoR):    
    data = pd.read_csv(archivoR,header=None, skiprows=0, index_col=0)#lee archivo X_rot.csv con las rotaciones
    arr = np.array(data)   #arreglo rotaciones por frame 
    f = arr[0]		    #fila 1 con nombres de marcadores
    return arr, f  #retorna arreglo de rotaciones y lista de headers de marcadores               


#funcion que recibe posiciones de la funcion offset, y rotaciones de funcion rotaciones
#retorna posiciones y rotaciones ajustados.
def ajustes(pos, rot):
    conta1, conta2, conta3, conta4 = 0, 0, 0, 0
    posiciones1, posiciones2, posiciones3 = [], [], []
    rotaciones = []
    v=int(len(rot[1])/3) #constante para el split 
    for i in rot[1:-1]:  #recorriendo rotaciones de A
        frame_ang = np.array_split(i,v)  #arreglo de angulos en grupos de 3
        frame_ang = np.array(frame_ang)  #se convirte a numpy
        Frame = frame_ang.astype(float)  # se convirte a floats con astype
        W = np.around(Frame, decimals=4) #se redondea
        rotaciones.insert(conta4, W)     #nueva lista rotaciones
        conta4 +=1
    Rot = rotaciones   #se guarda el arreglo Rot
    posi = list(pos)   #lista con posiciones
    for i in posi: #ahora se trabaja con offsets
        posiciones1.insert(conta1, i[1].strip()) #se inserta offsets en lista posiciones1
        conta1 +=1
    tiemp = np.array(posiciones1)  #se convierte la lista a numpy array
    for i in tiemp: 
        tmp = i.split(" ")         #se separa pos espacios como string
        posiciones2.insert(conta2, tmp)     #se inserta en lista posiciones2
        conta2 +=1 
    posiciones2 = np.array(posiciones2)     #se convierte a numpy   
    aux = posiciones2.astype(float)         #se convirte a float
    for i in aux:
        y = np.around(i, decimals=5)
        posiciones3.insert(conta3, y)       #se redondea y se inserta en lista posiciones3
        conta3 +=1
    posicionF = np.array_split(posiciones3,v)  #se separa en grupos de 3, en posicionF 
    return posicionF, Rot    #retorna rotaciones por frames, y posiciones

    
#funcion que recibe datos ajustados de angulos de Euler en secuencia ZXY y posicion XYZ
#retorna una posicion local
def matrices(angulos,pos):
    r = R.from_euler('zyx', angulos, degrees=True)  #secuencia ZXY y grados
    x = r.as_matrix() #obtiene matriz de rotacion de r
    y = r.apply(pos)  #multiplica matriz de rotacion y vector posicion.
    return y  #retorna posicion local


#recibe las posiciones de offset y rotaciones, ya ajustados
#retorna archivo localesX.csv con posiciones locales, ordenado en frames   
def sec(posi, rot, a, nombre1, nombre2):
    posFin = []
    cont1, cont2 = 0, 0
    v=int(len(rot[1])*3)      #constante para el reshape
    for j in rot:             #recorriendo arreglo de rotaciones, de la secuencia
        cont1 = 0
        for i in j:
            z = posi[cont1]   #recorriendo rotaciones, de cada frame 
            cont1 += 1
            w = matrices(i, z) #obteniendo posicion local, usando funcion matrices
            for h in w:
                posFin.insert(cont2, h) #creando arreglo con posiciones locales
                cont2 += 1
    cc = np.array(posFin)
    ccn= cc.reshape(-1,v) #ordena rotaciones
    kk= a.reshape(-1,v)   #ordena nombres de marcadores
    with open(nombre1, 'w', newline = '') as nombre2: #guardando en localesX.csv
        writer = csv.writer(nombre2)
        for i in kk:
            writer.writerow(i)  #escribe nombres 
        for i in ccn:
            writer.writerow(i)  #escribe posiciones locales
    return posFin 

def main():

    secuenciaA, secuenciaB = 'A.bvh', 'B.bvh'
    posA = offsets(secuenciaA)  #arreglo de offsets de A 
    posB = offsets(secuenciaB)  #arreglo de offsets de B 
   
    secRotA, secRotB = 'A_rot.csv', 'B_rot.csv'
    motionA, marcadoresA = rotaciones(secRotA)# 2 arreglos, rotaciones de A y nombres de marcadores
    motionB, marcadoresB = rotaciones(secRotB)# 2 arreglos, rotaciones de B y nombres de marcadores
    
    posiA, motionA = ajustes(posA, motionA)#arreglos de offsets y rotaciones de A ajustados.
    posiB, motionB = ajustes(posB, motionB)#arreglos de offsets y rotaciones de B ajustados.
    
    lA, lB, la, lb = 'localesA.csv','localesB.csv','localesA','localesB'
    posFinalA = sec(posiA, motionA, marcadoresA, lA, la) #Genera archivo CSV con posiciones locales de A
    posFinalB = sec(posiB, motionB, marcadoresB, lB, lb) #Genera archivo CSV con posiciones locales de B

main = main()


