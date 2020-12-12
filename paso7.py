#Step7.py "4D sequence time alignment" 
#Johnatan Quintero V. johnatan2425@hotmail.com

#Step7: Apartir de las 3 secuencias. Obtiene posiciones locales de cada una.
#       Obtiene el ISE para las curvas del promedio de orientacion, de cita y phi
#       Genera gráficas de error entre secuencias A y B. Y entre A y C.
#       Note que se utilizan funciones de pasos anteriores.

import math
import numpy as np    
import pandas as pd   
import csv              
import os                
import bvhtoolbox as bh          #libreria BVH
from scipy.spatial.transform import Rotation as R  #matrices de rotacion
from scipy import integrate      #area bajo la curva integral trapezoide
import matplotlib.pyplot as plt  #graficas

#genera el archivo C_rot.csv necesario para obtener posiciones locales
def rotacionesC():
    command = "bvh2csv -r ./C.bvh"   #se genera C_rot.csv
    print("Ejecutando el comando: \n", command)
    os.system(command)               #ejecuta comando en terminal  
              
           
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
    v=int(len(rot[1])*3)
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
    
     
#recibe archivo, de posiciones locales de una secuencia
#retorna posiciones locales en esfericas, cita y phi para cada marcador
def conversion(locales):
    pos = pd.read_csv(locales)  #lee archivo localesX.csv
    p = np.array(pos)
    esfer1, esfer2, esfer3 = [], [], []
    cont1, cont2, cont3 = 0, 0, 0
    for i in p:            #recorre cada frame, leyendo posiciones locales
        esfer1 = np.array_split(i, int(len(i)/3))
        cnt1=0
        for j in esfer1:
            cnt1+=3
            x=j[0]
            y=j[1]       #x,y,z de cada marcador
            z=j[2]
            aux1= math.pow(x, 2) + math.pow(y,2) #variable para ecuacion 1
            aux2= math.sqrt(aux1)  #ecuacion 1, variable para ecuacion 2
            if z!= 0: 
                aux3 = math.atan(aux2/z) #ecuacion 2,se obtiene cita en esfericas 
            else:
                z= z+0.00000001 #para no dividir entre cero
                aux3 = math.atan(aux2/z) #ecuacion 2,se obtiene cita en esfericas 
            z= z-0.00000001#se devuelve el valor
            if x != 0:
                aux4 = math.atan2(y,x) #ecse obtiene phi en esfericas
            else:
                x= x+0.00000001
                aux4 = math.atan2(y,x) #se obtiene componente phi en esfericas
            x= x-0.00000001
            cita= math.degrees(aux3) #pasando de radianes a grados
            phi= math.degrees(aux4)  
            esfer2.insert(cont1, round(cita, 5))
            cont1+=1
            esfer2.insert(cont1, round(phi, 5))
            cont1+=1
        cnt2=int((cnt1/3)*2)  #constante para el reshape
        esf2 = np.array(esfer2)
        es2 = esf2.reshape(-1, cnt2) #arreglo con posiciones locales en esfericas  
    return es2 #arreglo de posiciones locales en esfericas
    
    
#recibe arreglo de posiciones locales esfericas, y arreglos de tiempos de A, B y C
#devuelve area bajo la curva, de curva de latitud y longitud. Y dos arreglos de latitud y longitud. 
def validacion1(posEsf, tA, tB, tC):
    array1, array2 = [], []
    latitud, longitud = [], []
    cnt1, cnt2 = 0, 0
    tmp1= np.array([len(tA), len(tB), len(tC)])  #cantidad de frames de las 3 secuencias
    aux1= np.amin(tmp1)                          #guarda en aux1 la cantidad de frames, de la sec menor
    
    compare = np.array(posEsf[1:int(aux1-1)]) #arreglo con posiciones esfericas

    for i in compare: #recorriendo arreglo de posiciones locales en esfericas, cita y phi 
        array = np.array_split(i,int(len(i)/2))  #separando en grupos, para cada marcador
        lat1, lat2, lon1, lon2 = 0, 0, 0, 0
        for j in array[1:int(len(i)/2)]: #recorriendo marcadores
            lat1 = abs(j[0])     
            lon1 = abs(j[1])
            lat2 = lat2+lat1  #va sumando las latitudes de cada marcador
            lon2 = lon2+lon1  #va sumando las longitudes de cada marcador
        lat3= lat2/int(len(i)/2)  #promedio del total de suma de latitudes
        lon3= lon2/int(len(i)/2)  #promedio del total de suma de longitudes
        
        latitud.insert(cnt1, lat3)    #arreglo de promedio de latitudes
        longitud.insert(cnt1, lon3)   #arreglo de promedio de longitudes
        
        cnt1+=1
        
    reslat= (integrate.trapz(latitud, x=None))**2  #area bajo curva de latitud al cuadrado
    reslon= (integrate.trapz(longitud, x=None))**2 #area bajo curva de longitud al cuadrado
    return reslat, reslon, latitud, longitud    
 
 
#obtiene el ISE de curvas de latitud y longitud entre secuencias A y B. 
#genera graficas de curvas de latitud y longitud para validacion 
def validacionAB():
    Arot ,Brot, Crot = 'A_rot.csv','B_rot.csv','C_rot.csv'
    tA = leerTiempos(Arot)   #arreglos de tiempos de A
    tB = leerTiempos(Brot)   #arreglos de tiempos de B     
    tC = leerTiempos(Crot)   #arreglos de tiempos de C    
     
    localA, localB = 'localesA.csv','localesB.csv'   
    posEsfA = conversion(localA) #posiciones locales esfericas de A
    reslatA, reslonA, latitudA, longitudA = validacion1(posEsfA, tA, tB, tC)
    posEsfB = conversion(localB) #posiciones locales esfericas de B
    reslatB, reslonB, latitudB, longitudB = validacion1(posEsfB, tA, tB, tC)
    c=abs(reslatA-reslatB) #valor absoluto de la resta de las areas de curvas de latitud 
    cc=abs(reslonA-reslonB) #valor absoluto de la resta de las areas de curvas de longitud

    tmp1= np.array([len(tA), len(tB), len(tC)])
    aux1= np.amin(tmp1) #la cantidad de frames de la sec menor
    
    print('\nError entre secuencia A y secuencia B en los primeros',len(tA[0:aux1-2]),'frames')
    print('En curvas de latitud, ISE: ', int(c))
    print('En curvas de longitud, ISE: ', int(cc))
    
    
    t= tA[0:aux1-2] #vector de tiempo del tamano de la secuencia menor
    plt.plot(t, latitudA, 'r-',t, latitudB, 'b--') #grafica de latitudes de A y B
    plt.title('Curvas de latitud de secuencia A y B')
    plt.xlabel(r'Segundos, s')
    plt.ylabel(r'Latitud')
    plt.legend(('Secuencia A', 'Secuencia B'),prop = {'size': 10}, loc='upper right')
    plt.fill_between(t,latitudA, latitudB, color='grey')#rellenando area entre curvas
    plt.show()
   
    plt.plot(t, longitudA, 'r-',t, longitudB, 'b--') #grafica de longitudes de A y B
    plt.title('Curvas de longitud de secuencia A y B')
    plt.xlabel(r'Segundos, s')
    plt.ylabel(r'Longitud')
    plt.legend(('Secuencia A', 'Secuencia B'),prop = {'size': 10}, loc='upper right')
    plt.fill_between(t,longitudA, longitudB, color='grey')#rellenando area entre curvas
    plt.show()
    return int(c), int(cc)


#obtiene el ISE de curvas de latitud y longitud entre secuencias A y C. 
#genera graficas de curvas de latitud y longitud para validacion        
def validacionAC():
    Arot ,Brot, Crot = 'A_rot.csv','B_rot.csv','C_rot.csv'
    tA = leerTiempos(Arot)   #arreglos de tiempos de A
    tB = leerTiempos(Brot)   #arreglos de tiempos de B     
    tC = leerTiempos(Crot)   #arreglos de tiempos de C     
    
    localA, localC = 'localesA.csv','localesC.csv'   
    posEsfA = conversion(localA) #posiciones locales esfericas de A
    reslatA, reslonA, latitudA, longitudA = validacion1(posEsfA, tA, tB, tC)
    posEsfC = conversion(localC) #posiciones locales esfericas de C
    reslatC, reslonC, latitudC, longitudC = validacion1(posEsfC, tA, tB, tC)
    c=abs(reslatA-reslatC) #valor absoluto de la resta de las areas de curvas de latitud
    cc=abs(reslonA-reslonC) #valor absoluto de la resta de las areas de curvas de longitud

    tmp1= np.array([len(tA), len(tB), len(tC)])
    aux1= np.amin(tmp1) #la cantidad de frames de la sec menor
    
    print('\nError entre secuencia A y secuencia C en los primeros',len(tA[0:aux1-2]),'frames')
    print('En curvas de latitud, ISE: ', int(c))
    print('En curvas de longitud, ISE: ', int(cc))
    
    t= tA[0:aux1-2] #vector de tiempo del tamano de la secuencia menor
    plt.plot(t, latitudA, 'r-',t, latitudC, 'b--') #grafica de latitudes de A y C
    plt.title('Curvas de latitud de secuencia A y C')
    plt.xlabel(r'Segundos, s')
    plt.ylabel(r'Latitud')
    plt.legend(('Secuencia A', 'Secuencia C'),prop = {'size': 10}, loc='upper right')
    plt.fill_between(t,latitudA, latitudC, color='grey') #rellenando area entre curvas
    plt.show()
   
    plt.plot(t, longitudA, 'r-',t, longitudC, 'b--') #grafica de longitudes de A y C
    plt.title('Curvas de longitud de secuencia A y C')
    plt.xlabel(r'Segundos, s')
    plt.ylabel(r'Longitud')
    plt.legend(('Secuencia A', 'Secuencia C'),prop = {'size': 10}, loc='upper right')
    plt.fill_between(t,longitudA, longitudC, color='grey') #rellenando area entre curvas
    plt.show()
    return int(c), int(cc)


#recibe un archivo tipo X_rot.csv de una secuencia X
#devuelve un vector de tiempo de la secuencia    
def leerTiempos(archivo):
    tiemp = pd.read_csv(archivo)  
    t = np.array(tiemp)
    cont1 = 0
    tiempos = []
    for i in t:
        tiempos.insert(cont1, i[0])  #guardando los tiempos de frames 
        cont1+=1
    tiemp = np.array(tiempos)
    return tiemp  #arreglo de tiempos de frames


def main():
    rotacionesC()  #genera archivo C_rot.csv  
    secuenciaC = 'C.bvh'
    posC = offsets(secuenciaC)  #arreglo de offsets de C 
    
    secRotC = 'C_rot.csv'
    motionC, marcadoresC = rotaciones(secRotC)# 2 arreglos, rotaciones de C y nombres de marcadores
    posiC, motionC = ajustes(posC, motionC)#arreglos de offsets y rotaciones de C ajustados.
    
    lC, lc = 'localesC.csv','localesC'
    posFinalC = sec(posiC, motionC, marcadoresC, lC, lc) #Genera archivo CSV con posiciones locales de C

    cab1, cab2 = validacionAB() #ISE y grafica entre A y B
    cac1, cac2 = validacionAC() #ISE y grafica entre A y C
    ise1= abs(cab1-cac1)
    ise1= round((ise1*100)/cab1, 3)
    ise2= abs(cab2-cac2)
    ise2= round((ise2*100)/cab2, 3)
    print('Se mejoro un', ise1,'%','para latitud')
    print('Se mejoro un', ise2,'%','para longitud')
main = main()     


