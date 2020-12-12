#Step5.py "4D sequence time alignment" 
#Johnatan Quintero V. johnatan2425@hotmail.com

#Step5: Toma el string del step4, que tiene la secuencia alineada.
#       Donde los huecos son ceros. Y los numeros enteros son el frame respectivo.
#       Se alinean los frames de los archivos B_rot.csv y B_pos.csv   
#       Y se generan los archivos newB_rot.csv y newB_pos.csv. 
#       Necesarios para generar el BVH final.

from scipy.interpolate import lagrange #interpolacion de huecos
import numpy as np     
import pandas as pd    
import csv             
import math    

#recibe el nombre de un archivo X_rot.csv de rotaciones.
#devuelve un arreglo de los tiempos de cada frame
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


#recibe arreglo de tiempos de sec A y B, y lee el archivo textC.txt
#en el textC.txt hay 3 strings. Aca devuelve el string que tiene huecos
def alignment2(tiempoA, tiempoB):
    file = open('textC.txt', 'r')  #leyendo texto con alineacion
    sHuec1, sHuec2 = '', ''
    if len(tiempoA) > len(tiempoB):
        for line in file.readlines()[2]:    #si sec A es mayor, toma la sec B con huecos
            sHuec1 = line
            sHuec2 = sHuec2+sHuec1
    else:
        for line in file.readlines()[0]:   #si sec B es mayor, toma la sec A con huecos
            sHuec1 = line
            sHuec2 = sHuec2+sHuec1   
    file.close()
    return sHuec2    #devuelve un string de la sec con huecos 


#recibe un string con espacios como huecos, y los frames separados por una X
#devuelve string, con resultado con ceros como huecos y el numero de frame respectivo.    
def alignment3(huecos):
    cnt1, cnt2= 0, 0
    d1, d2= 0, 0
    arr= []
    for i in huecos:
        if i != 'X':        # X define un nuevo frame
            if i == '-':    #si hay un espacio suma 1 a cnt1
                cnt1+=1
            if cnt1 == 80:  #75 espacios que representan 1 hueco, 1 frame vacio
                d1+=1
                cnt1=0
        elif i == 'X':      # hay un nuevo frame
            cnt2+=1
            if d1 == 0:  #si no hay espacios
                arr.insert(d2,cnt2)  #inserta un numero que dice el frame respectivo
                d2+=1
            elif d1 != 0:#si hay espacios
                for i in range(d1):
                    arr.insert(d2, 0)#inserta un cero si hay minimo 75 espacios en blanco
                    d2+=1
                arr.insert(d2, cnt2) #inserta el frame que sigue, luego de poner los ceros 
                d2+=1      
            cnt1, d1 = 0, 0
    return arr      #devuelve arreglo con numeros para frames y ceros para huecos   


#recibe el string con ceros y numeros, recibe los vectores de tiempo de A y B
#devuelve newB_rot.csv y newB_pos.csv alineados. Si la secA es mayor que secB 
def alignmentAmayor(array, tA, tB, secB, new1, new2):  
    dB = pd.read_csv(secB,header=None, index_col=None)
    motionB = np.array(dB)     #arreglo con rotaciones del motion de B
    aux, cntt, cntf = 0, 0, 0
    head= list(motionB[0])     #cantidad de marcadores
    with open(new1, 'w', newline = '') as new2: #escribiendo en newB_rot.csv o en newB_pos.csv
        writer = csv.writer(new2)
        writer.writerow(motionB[0])  #escribe headers de marcadores
        for i in array:        #recorriendo string con numeros de frames y ceros 
            if i == 0 and cntf == 0:  #si hay un cero, y aun no llega el primer frame
                tmp1 = motionB[1]
                tmp1 = tmp1[1: len(head)]    #rellena antes del 1 frame
                tmp2 = np.array(round(tA[cntt],5))
                tmp3 = np.append(tmp2, tmp1)
                cntt+=1
                writer.writerow(tmp3)
            if i != 0:                #si es un numero
                tmp1 = motionB[i] #el frame que corresponce al entero que esta en i
                tmp1 = tmp1[1: len(head)]
                tmp2 = np.array(round(tA[cntt],5))
                tmp3 = np.append(tmp2, tmp1)
                writer.writerow(tmp3) #escribe el frame correspondiente
                cntf+=1 #aumenta contador de frames
                cntt+=1  #aumenta contador de tiempo
            if i == 0 and cntf!=0 and array[aux-1]!=0: #si es un cero y antes habia un numero          
                posi1, posi2 = [], []
                t1, t2, t3 = [], [], []
                c1, c2= 0, 0                
                posi1 = motionB[cntf]#frame actual
                posi1 = posi1[1: len(head)] #guarda pasado actual
                posi2 = motionB[cntf+1]#frame siguiente
                posi2 = posi2[1: len(head)] #guarda frame siguiente
                tiemp1 = round(tA[cntt-1],5) #tiempo pasado
                while array[aux+c1] == 0: #averiguando cuantos ceros vienen
                    c1+=1
                tiemp2 = round(tA[cntt+c1],5) #tiempo final
                while array[aux+c2] == 0:
                    c3=0
                    posi3 = []
                    for i in posi1:  #recorre marcadores del frame actual     
                        marcador1= round(float(i), 5)
                        marcador2= round(float(posi2[c3]), 5)
                        lagrangeX = np.append(tiemp1, tiemp2) #dos tiempos
                        lagrangeY = np.append(marcador1, marcador2) #dos puntos
                        interpol= lagrange(lagrangeX, lagrangeY) #interpolacion de lagrange
                        answ = interpol(round(tA[cntt], 5)) #evaluando la funcion en el tiempo que sigue
                        answ= round(answ, 5)
                        posi3.insert(c3, answ) #insertando en lista con resultados
                        c3+=1 
                    tp2 = np.array(round(tA[cntt],5)) #tiempo actual
                    tp3 = np.append(tp2, posi3) #frame para rellenar
                    c2+=1 #contador del while
                    cntt+=1 #aumenta contador tiempo
                    writer.writerow(tp3) #rellenando el renglon con el frame obtenido  
            aux+=1    


#recibe el string con ceros y numeros, recibe los vectores de tiempo de A y B
#devuelve newB_rot.csv y newB_pos.csv alineados. Si la secB es mayor que secA            
def alignmentBmayor(array, tA, tB, secB, new1, new2):
    dB = pd.read_csv(secB,header=None, index_col=None)
    motionB = np.array(dB)   #arreglo A de rotaciones por frame 
    aux1, aux2, aux3, cntf = 0, 0, 0, 0
    head= list(motionB[0])
    with open(new1, 'w', newline = '') as new2: #escribiendo en newB_rot.csv o en newB_pos.csv
        writer = csv.writer(new2)
        for i in array:  #recorriendo string con numeros de frames y ceros 
            if aux1 == 0:
                writer.writerow(motionB[0])  #escribe headers de marcadores
                aux1+=1
            if i != 0 and aux2 !=0:  #si es numero y es el primer frame
                tmp1 = motionB[aux2] #frame actual
                tmp1 = tmp1[1: len(head)]
                tmp2 = np.array(round(tA[aux3],5)) #tiempo actual
                aux3+=1
                tmp3 = np.append(tmp2, tmp1)
                writer.writerow(tmp3)  #escribe el frame actual
            if i != 0 and aux2 ==0:    #si es numero y ya no es el primer frame
                tmp1 = motionB[aux2+1] #frame siguiente
                tmp1 = tmp1[1: len(head)]
                tmp2 = np.array(round(tA[aux3],5)) #tiempo actual
                aux3+=1
                tmp3 = np.append(tmp2, tmp1)
                writer.writerow(tmp3)  #escribe el frame siguiente    
            aux2+=1 #note que si es cero no hace nada, porque se eliminan

def main():

    Arot,Brot, Bpos = 'A_rot.csv','B_rot.csv','B_pos.csv'
    tiempoA = leerTiempos(Arot)   #arreglos de tiempos de A
    tiempoB = leerTiempos(Brot)   #arreglos de tiempos de B
    
    secHuec = alignment2(tiempoA, tiempoB) #de los 3 strings de textC.txt, toma el que tiene espacios
    
    arr = alignment3(secHuec) #devuelve string de resultado con huecos como ceros y frames con numero
    
    dB = pd.read_csv('B_rot.csv',header=None, index_col=None)
    motionB = np.array(dB)     #arreglo con rotaciones del motion de B
    head= list(motionB[0])     
    
    newBrot, newBpos, newBr, newBp = 'newB_rot.csv','newB_pos.csv','newBrot','newBpos'
    if len(tiempoA) > len(tiempoB): #si A es mas largo que B
        alignmentAmayor(arr, tiempoA, tiempoB, Brot, newBrot, newBr) #genera newB_rot.csv
        alignmentAmayor(arr, tiempoA, tiempoB, Bpos, newBpos, newBp) #genera newB_pos.csv
    else:                           #si B es mas largo que A
        alignmentBmayor(arr, tiempoA, tiempoB, Brot, newBrot, newBr) #genera newB_rot.csv
        alignmentBmayor(arr, tiempoA, tiempoB, Bpos, newBpos, newBp) #genera newB_pos.csv

main =main()
