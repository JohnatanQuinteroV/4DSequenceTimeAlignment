#Step3.py "4D sequence time alignment" 
#Johnatan Quintero V. johnatan2425@hotmail.com

#Step3: Toma archivos localesA.csv y localesB.csv 
#       Convierte a cita y phi, en esfericas. Y realiza partición de esfera.
#       Se generan 2 archivos con string de codigos, textA.txt y textB.txt

import numpy as np    
import pandas as pd    
import csv             
import math      
       
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
                aux4 = math.atan2(y,x) #se obtiene phi en esfericas
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


#recibe posiciones locales esfericas, de la funcion conversión.    
#retorna string de codigo para la secuencia. Usando particion de esfera, en áreas iguales. 
def particion(es2):
    array, arrayAux = [], []
    codigo = []
    cnt1, cnt2 = 0, 0
    esfera = np.array(es2)
    for i in esfera: #recorriendo arreglo de posiciones locales en esfericas, cita y phi
        array = np.array_split(i,int(len(i)/2))  #separando en grupos, para cada marcador
        for j in array[1:int(len(i)/2)]: #recorriendo arreglo de cada frame
            if 80.1375 < j[0] <= 90: #secciones de latitudes
                codigo.insert(cnt1,'AA')
                cnt1+=1
                if 0 <= j[1] <= 120: #secciones de longitud, segun la latitud
                    codigo.insert(cnt1,'ba')
                    cnt1+=1
                elif 120 < j[1] <= 180:
                    codigo.insert(cnt1,'bb')
                    cnt1+=1
                elif -180 <= j[1] <= -120:
                    codigo.insert(cnt1,'bb')
                    cnt1+=1
                elif -120 < j[1] < 0:
                    codigo.insert(cnt1,'bc')
                    cnt1+=1
            elif 70.2010 < j[0] <= 80.1375:  #latitud
                codigo.insert(cnt1,'AB')
                cnt1+=1
                if 0 <= j[1] <= 40:#seccion de longitudes
                    codigo.insert(cnt1,'ba')
                    cnt1+=1
                elif 40 < j[1] <= 80:
                    codigo.insert(cnt1,'bb')
                    cnt1+=1
                elif 80 < j[1] <= 120:
                    codigo.insert(cnt1,'bc')
                    cnt1+=1
                elif 120 < j[1] <= 160:
                    codigo.insert(cnt1,'bd')
                    cnt1+=1
                elif 160 < j[1] <= 180:
                    codigo.insert(cnt1,'be')
                    cnt1+=1
                elif -180 <= j[1] <= -160:
                    codigo.insert(cnt1,'bf')
                    cnt1+=1
                elif -160 < j[1] <= -120:
                    codigo.insert(cnt1,'bg')
                    cnt1+=1
                elif -120 < j[1] <= -80:
                    codigo.insert(cnt1,'bh')
                    cnt1+=1
                elif -80 < j[1] <= -40:
                    codigo.insert(cnt1,'bi')
                    cnt1+=1
                elif -40 < j[1] < 0:
                    codigo.insert(cnt1,'bj')
                    cnt1+=1
                    
            elif 60.1113 < j[0] <= 70.2010:#latitud
                codigo.insert(cnt1,'AC')
                cnt1+=1
                if 0 <= j[1] <= 24:#longitudes
                    codigo.insert(cnt1,'ba')
                    cnt1+=1
                elif 24 < j[1] <= 48:
                    codigo.insert(cnt1,'bb')
                    cnt1+=1
                elif 48 < j[1] <= 72:
                    codigo.insert(cnt1,'bc')
                    cnt1+=1
                elif 72 < j[1] <= 96:
                    codigo.insert(cnt1,'bd')
                    cnt1+=1
                elif 96 < j[1] <= 120:
                    codigo.insert(cnt1,'be')
                    cnt1+=1
                elif 120 < j[1] <= 144:
                    codigo.insert(cnt1,'bf')
                    cnt1+=1
                elif 144 < j[1] <= 168:
                    codigo.insert(cnt1,'bg')
                    cnt1+=1
                elif 168 < j[1] <= 180:
                    codigo.insert(cnt1,'bh')
                    cnt1+=1
                elif -180 <= j[1] <= -168:
                    codigo.insert(cnt1,'bh')
                    cnt1+=1
                elif -168 < j[1] <= -144:
                    codigo.insert(cnt1,'bi')
                    cnt1+=1
                elif -144 < j[1] <= -120:
                    codigo.insert(cnt1,'bj')
                    cnt1+=1
                elif -120 < j[1] <= -96:
                    codigo.insert(cnt1,'bk')
                    cnt1+=1
                elif -96 < j[1] <= -72:
                    codigo.insert(cnt1,'bl')
                    cnt1+=1
                elif -72 < j[1] <= -48:
                    codigo.insert(cnt1,'bm')
                    cnt1+=1
                elif -48 < j[1] <= -24:
                    codigo.insert(cnt1,'bn')
                    cnt1+=1
                elif -24 < j[1] < 0:
                    codigo.insert(cnt1,'bo')
                    cnt1+=1
            elif 50.2170 < j[0] <= 60.1113:
                codigo.insert(cnt1,'AD')
                cnt1+=1
                if 0 <= j[1] <= 18:
                    codigo.insert(cnt1,'ba')
                    cnt1+=1
                elif 18 < j[1] <= 36:
                    codigo.insert(cnt1,'bb')
                    cnt1+=1
                elif 36 < j[1] <= 54:
                    codigo.insert(cnt1,'bc')
                    cnt1+=1
                elif 54 < j[1] <= 72:
                    codigo.insert(cnt1,'bd')
                    cnt1+=1
                elif 72 < j[1] <= 90:
                    codigo.insert(cnt1,'be')
                    cnt1+=1
                elif 90 < j[1] <= 108:
                    codigo.insert(cnt1,'bf')
                    cnt1+=1
                elif 108 < j[1] <= 126:
                    codigo.insert(cnt1,'bg')
                    cnt1+=1
                elif 126 < j[1] <= 144:
                    codigo.insert(cnt1,'bh')
                    cnt1+=1
                elif 144 < j[1] <= 162:
                    codigo.insert(cnt1,'bi')
                    cnt1+=1
                elif 162 < j[1] <= 180:
                    codigo.insert(cnt1,'bj')
                    cnt1+=1
                elif -180 <= j[1] <= -162:
                    codigo.insert(cnt1,'bk')
                    cnt1+=1
                elif -162 < j[1] <= -144:
                    codigo.insert(cnt1,'bl')
                    cnt1+=1
                elif -144 < j[1] <= -126:
                    codigo.insert(cnt1,'bm')
                    cnt1+=1
                elif -126 < j[1] <= -108:
                    codigo.insert(cnt1,'bn')
                    cnt1+=1
                elif -108 < j[1] <= -90:
                    codigo.insert(cnt1,'bo')
                    cnt1+=1
                elif -90 < j[1] <= -72:
                    codigo.insert(cnt1,'bp')
                    cnt1+=1
                elif -72 < j[1] <= -54:
                    codigo.insert(cnt1,'bq')
                    cnt1+=1
                elif -54 < j[1] <= -36:
                    codigo.insert(cnt1,'br')
                    cnt1+=1
                elif -36 < j[1] <= -18:
                    codigo.insert(cnt1,'bs')
                    cnt1+=1
                elif -18 < j[1] < 0:
                    codigo.insert(cnt1,'bt')
                    cnt1+=1
            elif 40.5602 < j[0] <= 50.2170:
                codigo.insert(cnt1,'AE')
                cnt1+=1
                if 0 <= j[1] <= 15:
                    codigo.insert(cnt1,'ba')
                    cnt1+=1
                elif 15 < j[1] <= 30:
                    codigo.insert(cnt1,'bb')
                    cnt1+=1
                elif 30 < j[1] <= 45:
                    codigo.insert(cnt1,'bc')
                    cnt1+=1
                elif 45 < j[1] <= 60:
                    codigo.insert(cnt1,'bd')
                    cnt1+=1
                elif 60 < j[1] <= 75:
                    codigo.insert(cnt1,'be')
                    cnt1+=1
                elif 75 < j[1] <= 90:
                    codigo.insert(cnt1,'bf')
                    cnt1+=1
                elif 90 < j[1] <= 105:
                    codigo.insert(cnt1,'bg')
                    cnt1+=1
                elif 105 < j[1] <= 120:
                    codigo.insert(cnt1,'bh')
                    cnt1+=1
                elif 120 < j[1] <= 135:
                    codigo.insert(cnt1,'bi')
                    cnt1+=1
                elif 135 < j[1] <= 150:
                    codigo.insert(cnt1,'bj')
                    cnt1+=1
                elif 150 < j[1] <= 165:
                    codigo.insert(cnt1,'bk')
                    cnt1+=1
                elif 165 < j[1] <= 180:
                    codigo.insert(cnt1,'bl')
                    cnt1+=1
                elif -180 <= j[1] <= -165:
                    codigo.insert(cnt1,'bm')
                    cnt1+=1
                elif -165 < j[1] <= -150:
                    codigo.insert(cnt1,'bn')
                    cnt1+=1
                elif -150 < j[1] <= -135:
                    codigo.insert(cnt1,'bo')
                    cnt1+=1
                elif -135 < j[1] <= -120:
                    codigo.insert(cnt1,'bp')
                    cnt1+=1
                elif -120 < j[1] <= -105:
                    codigo.insert(cnt1,'bq')
                    cnt1+=1
                elif -105 < j[1] <= -90:
                    codigo.insert(cnt1,'br')
                    cnt1+=1
                elif -90 < j[1] <= -75:
                    codigo.insert(cnt1,'bs')
                    cnt1+=1
                elif -75 < j[1] <= -60:
                    codigo.insert(cnt1,'bt')
                    cnt1+=1
                elif -60 < j[1] <= -45:
                    codigo.insert(cnt1,'bu')
                    cnt1+=1
                elif -45 < j[1] <= -30:
                    codigo.insert(cnt1,'bv')
                    cnt1+=1
                elif -30 < j[1] <= -15:
                    codigo.insert(cnt1,'bw')
                    cnt1+=1
                elif -15 < j[1] < 0:
                    codigo.insert(cnt1,'bx')
                    cnt1+=1
            elif 30.1631 < j[0] <= 40.5602:
                codigo.insert(cnt1,'AF')
                cnt1+=1                
                if 0 <= j[1] <= 12:
                    codigo.insert(cnt1,'ba')
                    cnt1+=1
                elif 12 < j[1] <= 24:
                    codigo.insert(cnt1,'bb')
                    cnt1+=1
                elif 24 < j[1] <= 36:
                    codigo.insert(cnt1,'bc')
                    cnt1+=1
                elif 36 < j[1] <= 48:
                    codigo.insert(cnt1,'bd')
                    cnt1+=1
                elif 48 < j[1] <= 60:
                    codigo.insert(cnt1,'be')
                    cnt1+=1
                elif 60 < j[1] <= 72:
                    codigo.insert(cnt1,'bf')
                    cnt1+=1
                elif 72 < j[1] <= 84:
                    codigo.insert(cnt1,'bg')
                    cnt1+=1
                elif 84 < j[1] <= 96:
                    codigo.insert(cnt1,'bh')
                    cnt1+=1
                elif 96 < j[1] <= 108:
                    codigo.insert(cnt1,'bi')
                    cnt1+=1
                elif 108 < j[1] <= 120:
                    codigo.insert(cnt1,'bj')
                    cnt1+=1
                elif 120 < j[1] <= 132:
                    codigo.insert(cnt1,'bk')
                    cnt1+=1
                elif 132 < j[1] <= 144:
                    codigo.insert(cnt1,'bl')
                    cnt1+=1
                elif 144 < j[1] <= 156:
                    codigo.insert(cnt1,'bm')
                    cnt1+=1
                elif 156 < j[1] <= 168:
                    codigo.insert(cnt1,'bn')
                    cnt1+=1
                elif 168 < j[1] <= 180:
                    codigo.insert(cnt1,'bo')
                    cnt1+=1
                elif -180 <= j[1] <= -168:
                    codigo.insert(cnt1,'bp')
                    cnt1+=1
                elif -168 < j[1] <= -156:
                    codigo.insert(cnt1,'bq')
                    cnt1+=1
                elif -156 < j[1] <= -144:
                    codigo.insert(cnt1,'br')
                    cnt1+=1
                elif -144 < j[1] <= -132:
                    codigo.insert(cnt1,'bs')
                    cnt1+=1
                elif -132 < j[1] <= -120:
                    codigo.insert(cnt1,'bt')
                    cnt1+=1
                elif -120 < j[1] <= -108:
                    codigo.insert(cnt1,'bu')
                    cnt1+=1
                elif -108 < j[1] <= -96:
                    codigo.insert(cnt1,'bv')
                    cnt1+=1
                elif -96 < j[1] <= -84:
                    codigo.insert(cnt1,'bw')
                    cnt1+=1
                elif -84 < j[1] <= -72:
                    codigo.insert(cnt1,'bx')
                    cnt1+=1
                elif -72 < j[1] <= 60:
                    codigo.insert(cnt1,'by')
                    cnt1+=1
                elif -60 < j[1] <= -48:
                    codigo.insert(cnt1,'bz')
                    cnt1+=1
                elif -48 < j[1] <= -36:
                    codigo.insert(cnt1,'zy')
                    cnt1+=1
                elif -36 < j[1] <= -24:
                    codigo.insert(cnt1,'xw')
                    cnt1+=1
                elif -24 < j[1] <= -12:
                    codigo.insert(cnt1,'vu')
                    cnt1+=1
                elif -12 < j[1] < 0:
                    codigo.insert(cnt1,'ts')
                    cnt1+=1
            elif 20.7738 < j[0] <= 30.1631:
                codigo.insert(cnt1,'AG')
                cnt1+=1                 
                if 0 <= j[1] <= 12:
                    codigo.insert(cnt1,'ba')
                    cnt1+=1
                elif 12 < j[1] <= 24:
                    codigo.insert(cnt1,'bb')
                    cnt1+=1
                elif 24 < j[1] <= 36:
                    codigo.insert(cnt1,'bc')
                    cnt1+=1
                elif 36 < j[1] <= 48:
                    codigo.insert(cnt1,'bd')
                    cnt1+=1
                elif 48 < j[1] <= 60:
                    codigo.insert(cnt1,'be')
                    cnt1+=1
                elif 60 < j[1] <= 72:
                    codigo.insert(cnt1,'bf')
                    cnt1+=1
                elif 72 < j[1] <= 84:
                    codigo.insert(cnt1,'bg')
                    cnt1+=1
                elif 84 < j[1] <= 96:
                    codigo.insert(cnt1,'bh')
                    cnt1+=1
                elif 96 < j[1] <= 108:
                    codigo.insert(cnt1,'bi')
                    cnt1+=1
                elif 108 < j[1] <= 120:
                    codigo.insert(cnt1,'bj')
                    cnt1+=1
                elif 120 < j[1] <= 132:
                    codigo.insert(cnt1,'bk')
                    cnt1+=1
                elif 132 < j[1] <= 144:
                    codigo.insert(cnt1,'bl')
                    cnt1+=1
                elif 144 < j[1] <= 156:
                    codigo.insert(cnt1,'bm')
                    cnt1+=1
                elif 156 < j[1] <= 168:
                    codigo.insert(cnt1,'bn')
                    cnt1+=1
                elif 168 < j[1] <= 180:
                    codigo.insert(cnt1,'bo')
                    cnt1+=1
                elif -180 <= j[1] <= -168:
                    codigo.insert(cnt1,'bp')
                    cnt1+=1
                elif -168 < j[1] <= -156:
                    codigo.insert(cnt1,'bq')
                    cnt1+=1
                elif -156 < j[1] <= -144:
                    codigo.insert(cnt1,'br')
                    cnt1+=1
                elif -144 < j[1] <= -132:
                    codigo.insert(cnt1,'bs')
                    cnt1+=1
                elif -132 < j[1] <= -120:
                    codigo.insert(cnt1,'bt')
                    cnt1+=1
                elif -120 < j[1] <= -108:
                    codigo.insert(cnt1,'bu')
                    cnt1+=1
                elif -108 < j[1] <= -96:
                    codigo.insert(cnt1,'bv')
                    cnt1+=1
                elif -96 < j[1] <= -84:
                    codigo.insert(cnt1,'bw')
                    cnt1+=1
                elif -84 < j[1] <= -72:
                    codigo.insert(cnt1,'bx')
                    cnt1+=1
                elif -72 < j[1] <= 60:
                    codigo.insert(cnt1,'by')
                    cnt1+=1
                elif -60 < j[1] <= -48:
                    codigo.insert(cnt1,'bz')
                    cnt1+=1
                elif -48 < j[1] <= -36:
                    codigo.insert(cnt1,'zy')
                    cnt1+=1
                elif -36 < j[1] <= -24:
                    codigo.insert(cnt1,'xw')
                    cnt1+=1
                elif -24 < j[1] <= -12:
                    codigo.insert(cnt1,'vu')
                    cnt1+=1
                elif -12 < j[1] < 0:
                    codigo.insert(cnt1,'ts')
                    cnt1+=1
            elif 10.2148 < j[0] <= 20.7738:
                codigo.insert(cnt1,'AH')
                cnt1+=1                
                if 0 <= j[1] <= 10:
                    codigo.insert(cnt1,'ba')
                    cnt1+=1
                elif 10 < j[1] <= 20:
                    codigo.insert(cnt1,'bb')
                    cnt1+=1
                elif 20 < j[1] <= 30:
                    codigo.insert(cnt1,'bc')
                    cnt1+=1
                elif 30 < j[1] <= 40:
                    codigo.insert(cnt1,'bd')
                    cnt1+=1
                elif 40 < j[1] <= 50:
                    codigo.insert(cnt1,'be')
                    cnt1+=1
                elif 50 < j[1] <= 60:
                    codigo.insert(cnt1,'bf')
                    cnt1+=1
                elif 60 < j[1] <= 70:
                    codigo.insert(cnt1,'bg')
                    cnt1+=1
                elif 70 < j[1] <= 80:
                    codigo.insert(cnt1,'bh')
                    cnt1+=1
                elif 80 < j[1] <= 90:
                    codigo.insert(cnt1,'bi')
                    cnt1+=1
                elif 90 < j[1] <= 100:
                    codigo.insert(cnt1,'bj')
                    cnt1+=1
                elif 100 < j[1] <= 110:
                    codigo.insert(cnt1,'bk')
                    cnt1+=1
                elif 110 < j[1] <= 120:
                    codigo.insert(cnt1,'bl')
                    cnt1+=1
                elif 120 < j[1] <= 130:
                    codigo.insert(cnt1,'bm')
                    cnt1+=1
                elif 130 < j[1] <= 140:
                    codigo.insert(cnt1,'bn')
                    cnt1+=1
                elif 140 < j[1] <= 150:
                    codigo.insert(cnt1,'bo')
                    cnt1+=1
                elif 150 < j[1] <= 160:
                    codigo.insert(cnt1,'bp')
                    cnt1+=1
                elif 160 < j[1] <= 170:
                    codigo.insert(cnt1,'bq')
                    cnt1+=1
                elif 170 < j[1] <= 180:
                    codigo.insert(cnt1,'br')
                    cnt1+=1
                elif -180 <= j[1] <= -170:
                    codigo.insert(cnt1,'bs')
                    cnt1+=1
                elif -170 < j[1] <= -160:
                    codigo.insert(cnt1,'bt')
                    cnt1+=1
                elif -150 < j[1] <= -140:
                    codigo.insert(cnt1,'bu')
                    cnt1+=1
                elif -140 < j[1] <= -130:
                    codigo.insert(cnt1,'bv')
                    cnt1+=1
                elif -130 < j[1] <= -120:
                    codigo.insert(cnt1,'bw')
                    cnt1+=1
                elif -120 < j[1] <= -110:
                    codigo.insert(cnt1,'bx')
                    cnt1+=1
                elif -110 < j[1] <= -100:
                    codigo.insert(cnt1,'by')
                    cnt1+=1
                elif -100 < j[1] <= -90:
                    codigo.insert(cnt1,'bz')
                    cnt1+=1
                elif -90 < j[1] <= -80:
                    codigo.insert(cnt1,'zy')
                    cnt1+=1
                elif -80 < j[1] <= -70:
                    codigo.insert(cnt1,'xw')
                    cnt1+=1
                elif -70 < j[1] <= -60:
                    codigo.insert(cnt1,'vu')
                    cnt1+=1
                elif -60 < j[1] <= -50:
                    codigo.insert(cnt1,'ts')
                    cnt1+=1
                elif -50 < j[1] <= -40:
                    codigo.insert(cnt1,'rq')
                    cnt1+=1
                elif -40 < j[1] <= -30:
                    codigo.insert(cnt1,'po')
                    cnt1+=1
                elif -30 < j[1] <= -20:
                    codigo.insert(cnt1,'nm')
                    cnt1+=1
                elif -20 < j[1] <= -10:
                    codigo.insert(cnt1,'lk')
                    cnt1+=1
                elif -10 < j[1] < 0:
                    codigo.insert(cnt1,'ji')
                    cnt1+=1

            elif 0 <= j[0] <= 10.2148:
                codigo.insert(cnt1,'AI')
                cnt1+=1
                if 0 <= j[1] <= 10:
                    codigo.insert(cnt1,'ba')
                    cnt1+=1
                elif 10 < j[1] <= 20:
                    codigo.insert(cnt1,'bb')
                    cnt1+=1
                elif 20 < j[1] <= 30:
                    codigo.insert(cnt1,'bc')
                    cnt1+=1
                elif 30 < j[1] <= 40:
                    codigo.insert(cnt1,'bd')
                    cnt1+=1
                elif 40 < j[1] <= 50:
                    codigo.insert(cnt1,'be')
                    cnt1+=1
                elif 50 < j[1] <= 60:
                    codigo.insert(cnt1,'bf')
                    cnt1+=1
                elif 60 < j[1] <= 70:
                    codigo.insert(cnt1,'bg')
                    cnt1+=1
                elif 70 < j[1] <= 80:
                    codigo.insert(cnt1,'bh')
                    cnt1+=1
                elif 80 < j[1] <= 90:
                    codigo.insert(cnt1,'bi')
                    cnt1+=1
                elif 90 < j[1] <= 100:
                    codigo.insert(cnt1,'bj')
                    cnt1+=1
                elif 100 < j[1] <= 110:
                    codigo.insert(cnt1,'bk')
                    cnt1+=1
                elif 110 < j[1] <= 120:
                    codigo.insert(cnt1,'bl')
                    cnt1+=1
                elif 120 < j[1] <= 130:
                    codigo.insert(cnt1,'bm')
                    cnt1+=1
                elif 130 < j[1] <= 140:
                    codigo.insert(cnt1,'bn')
                    cnt1+=1
                elif 140 < j[1] <= 150:
                    codigo.insert(cnt1,'bo')
                    cnt1+=1
                elif 150 < j[1] <= 160:
                    codigo.insert(cnt1,'bp')
                    cnt1+=1
                elif 160 < j[1] <= 170:
                    codigo.insert(cnt1,'bq')
                    cnt1+=1
                elif 170 < j[1] <= 180:
                    codigo.insert(cnt1,'br')
                    cnt1+=1
                elif -180 <= j[1] <= -170:
                    codigo.insert(cnt1,'bs')
                    cnt1+=1
                elif -170 < j[1] <= -160:
                    codigo.insert(cnt1,'bt')
                    cnt1+=1
                elif -150 < j[1] <= -140:
                    codigo.insert(cnt1,'bu')
                    cnt1+=1
                elif -140 < j[1] <= -130:
                    codigo.insert(cnt1,'bv')
                    cnt1+=1
                elif -130 < j[1] <= -120:
                    codigo.insert(cnt1,'bw')
                    cnt1+=1
                elif -120 < j[1] <= -110:
                    codigo.insert(cnt1,'bx')
                    cnt1+=1
                elif -110 < j[1] <= -100:
                    codigo.insert(cnt1,'by')
                    cnt1+=1
                elif -100 < j[1] <= -90:
                    codigo.insert(cnt1,'bz')
                    cnt1+=1
                elif -90 < j[1] <= -80:
                    codigo.insert(cnt1,'zy')
                    cnt1+=1
                elif -80 < j[1] <= -70:
                    codigo.insert(cnt1,'xw')
                    cnt1+=1
                elif -70 < j[1] <= -60:
                    codigo.insert(cnt1,'vu')
                    cnt1+=1
                elif -60 < j[1] <= -50:
                    codigo.insert(cnt1,'ts')
                    cnt1+=1
                elif -50 < j[1] <= -40:
                    codigo.insert(cnt1,'rq')
                    cnt1+=1
                elif -40 < j[1] <= -30:
                    codigo.insert(cnt1,'po')
                    cnt1+=1
                elif -30 < j[1] <= -20:
                    codigo.insert(cnt1,'nm')
                    cnt1+=1
                elif -20 < j[1] <= -10:
                    codigo.insert(cnt1,'lk')
                    cnt1+=1
                elif -10 < j[1] < 0:
                    codigo.insert(cnt1,'ji')
                    cnt1+=1
                    
            if -80.1375 > j[0] >= -90:
                codigo.insert(cnt1,'AJ')
                cnt1+=1
                if 0 <= j[1] <= 120:
                    codigo.insert(cnt1,'ba')
                    cnt1+=1
                elif 120 < j[1] <= 180:
                    codigo.insert(cnt1,'bb')
                    cnt1+=1
                elif -180 <= j[1] <= -120:
                    codigo.insert(cnt1,'bb')
                    cnt1+=1
                elif -120 < j[1] < 0:
                    codigo.insert(cnt1,'bc')
                    cnt1+=1
            elif -70.2010 > j[0] >= -80.1375:
                codigo.insert(cnt1,'AK')
                cnt1+=1
                if 0 <= j[1] <= 40:
                    codigo.insert(cnt1,'ba')
                    cnt1+=1
                elif 40 < j[1] <= 80:
                    codigo.insert(cnt1,'bb')
                    cnt1+=1
                elif 80 < j[1] <= 120:
                    codigo.insert(cnt1,'bc')
                    cnt1+=1
                elif 120 < j[1] <= 160:
                    codigo.insert(cnt1,'bd')
                    cnt1+=1
                elif 160 < j[1] <= 180:
                    codigo.insert(cnt1,'be')
                    cnt1+=1
                elif -180 <= j[1] <= -160:
                    codigo.insert(cnt1,'bf')
                    cnt1+=1
                elif -160 < j[1] <= -120:
                    codigo.insert(cnt1,'bg')
                    cnt1+=1
                elif -120 < j[1] <= -80:
                    codigo.insert(cnt1,'bh')
                    cnt1+=1
                elif -80 < j[1] <= -40:
                    codigo.insert(cnt1,'bi')
                    cnt1+=1
                elif -40 < j[1] < 0:
                    codigo.insert(cnt1,'bj')
                    cnt1+=1
                    
            elif -60.1113 > j[0] >= -70.2010:
                codigo.insert(cnt1,'AL')
                cnt1+=1
                if 0 <= j[1] <= 24:
                    codigo.insert(cnt1,'ba')
                    cnt1+=1
                elif 24 < j[1] <= 48:
                    codigo.insert(cnt1,'bb')
                    cnt1+=1
                elif 48 < j[1] <= 72:
                    codigo.insert(cnt1,'bc')
                    cnt1+=1
                elif 72 < j[1] <= 96:
                    codigo.insert(cnt1,'bd')
                    cnt1+=1
                elif 96 < j[1] <= 120:
                    codigo.insert(cnt1,'be')
                    cnt1+=1
                elif 120 < j[1] <= 144:
                    codigo.insert(cnt1,'bf')
                    cnt1+=1
                elif 144 < j[1] <= 168:
                    codigo.insert(cnt1,'bg')
                    cnt1+=1
                elif 168 < j[1] <= 180:
                    codigo.insert(cnt1,'bh')
                    cnt1+=1
                elif -180 <= j[1] <= -168:
                    codigo.insert(cnt1,'bh')
                    cnt1+=1
                elif -168 < j[1] <= -144:
                    codigo.insert(cnt1,'bi')
                    cnt1+=1
                elif -144 < j[1] <= -120:
                    codigo.insert(cnt1,'bj')
                    cnt1+=1
                elif -120 < j[1] <= -96:
                    codigo.insert(cnt1,'bk')
                    cnt1+=1
                elif -96 < j[1] <= -72:
                    codigo.insert(cnt1,'bl')
                    cnt1+=1
                elif -72 < j[1] <= -48:
                    codigo.insert(cnt1,'bm')
                    cnt1+=1
                elif -48 < j[1] <= -24:
                    codigo.insert(cnt1,'bn')
                    cnt1+=1
                elif -24 < j[1] < 0:
                    codigo.insert(cnt1,'bo')
                    cnt1+=1
            elif -50.2170 > j[0] >= -60.1113:
                codigo.insert(cnt1,'AM')
                cnt1+=1
                if 0 <= j[1] <= 18:
                    codigo.insert(cnt1,'ba')
                    cnt1+=1
                elif 18 < j[1] <= 36:
                    codigo.insert(cnt1,'bb')
                    cnt1+=1
                elif 36 < j[1] <= 54:
                    codigo.insert(cnt1,'bc')
                    cnt1+=1
                elif 54 < j[1] <= 72:
                    codigo.insert(cnt1,'bd')
                    cnt1+=1
                elif 72 < j[1] <= 90:
                    codigo.insert(cnt1,'be')
                    cnt1+=1
                elif 90 < j[1] <= 108:
                    codigo.insert(cnt1,'bf')
                    cnt1+=1
                elif 108 < j[1] <= 126:
                    codigo.insert(cnt1,'bg')
                    cnt1+=1
                elif 126 < j[1] <= 144:
                    codigo.insert(cnt1,'bh')
                    cnt1+=1
                elif 144 < j[1] <= 162:
                    codigo.insert(cnt1,'bi')
                    cnt1+=1
                elif 162 < j[1] <= 180:
                    codigo.insert(cnt1,'bj')
                    cnt1+=1
                elif -180 <= j[1] <= -162:
                    codigo.insert(cnt1,'bk')
                    cnt1+=1
                elif -162 < j[1] <= -144:
                    codigo.insert(cnt1,'bl')
                    cnt1+=1
                elif -144 < j[1] <= -126:
                    codigo.insert(cnt1,'bm')
                    cnt1+=1
                elif -126 < j[1] <= -108:
                    codigo.insert(cnt1,'bn')
                    cnt1+=1
                elif -108 < j[1] <= -90:
                    codigo.insert(cnt1,'bo')
                    cnt1+=1
                elif -90 < j[1] <= -72:
                    codigo.insert(cnt1,'bp')
                    cnt1+=1
                elif -72 < j[1] <= -54:
                    codigo.insert(cnt1,'bq')
                    cnt1+=1
                elif -54 < j[1] <= -36:
                    codigo.insert(cnt1,'br')
                    cnt1+=1
                elif -36 < j[1] <= -18:
                    codigo.insert(cnt1,'bs')
                    cnt1+=1
                elif -18 < j[1] < 0:
                    codigo.insert(cnt1,'bt')
                    cnt1+=1
            elif -40.5602 > j[0] >= -50.2170:
                codigo.insert(cnt1,'AN')
                cnt1+=1
                if 0 <= j[1] <= 15:
                    codigo.insert(cnt1,'ba')
                    cnt1+=1
                elif 15 < j[1] <= 30:
                    codigo.insert(cnt1,'bb')
                    cnt1+=1
                elif 30 < j[1] <= 45:
                    codigo.insert(cnt1,'bc')
                    cnt1+=1
                elif 45 < j[1] <= 60:
                    codigo.insert(cnt1,'bd')
                    cnt1+=1
                elif 60 < j[1] <= 75:
                    codigo.insert(cnt1,'be')
                    cnt1+=1
                elif 75 < j[1] <= 90:
                    codigo.insert(cnt1,'bf')
                    cnt1+=1
                elif 90 < j[1] <= 105:
                    codigo.insert(cnt1,'bg')
                    cnt1+=1
                elif 105 < j[1] <= 120:
                    codigo.insert(cnt1,'bh')
                    cnt1+=1
                elif 120 < j[1] <= 135:
                    codigo.insert(cnt1,'bi')
                    cnt1+=1
                elif 135 < j[1] <= 150:
                    codigo.insert(cnt1,'bj')
                    cnt1+=1
                elif 150 < j[1] <= 165:
                    codigo.insert(cnt1,'bk')
                    cnt1+=1
                elif 165 < j[1] <= 180:
                    codigo.insert(cnt1,'bl')
                    cnt1+=1
                elif -180 <= j[1] <= -165:
                    codigo.insert(cnt1,'bm')
                    cnt1+=1
                elif -165 < j[1] <= -150:
                    codigo.insert(cnt1,'bn')
                    cnt1+=1
                elif -150 < j[1] <= -135:
                    codigo.insert(cnt1,'bo')
                    cnt1+=1
                elif -135 < j[1] <= -120:
                    codigo.insert(cnt1,'bp')
                    cnt1+=1
                elif -120 < j[1] <= -105:
                    codigo.insert(cnt1,'bq')
                    cnt1+=1
                elif -105 < j[1] <= -90:
                    codigo.insert(cnt1,'br')
                    cnt1+=1
                elif -90 < j[1] <= -75:
                    codigo.insert(cnt1,'bs')
                    cnt1+=1
                elif -75 < j[1] <= -60:
                    codigo.insert(cnt1,'bt')
                    cnt1+=1
                elif -60 < j[1] <= -45:
                    codigo.insert(cnt1,'bu')
                    cnt1+=1
                elif -45 < j[1] <= -30:
                    codigo.insert(cnt1,'bv')
                    cnt1+=1
                elif -30 < j[1] <= -15:
                    codigo.insert(cnt1,'bw')
                    cnt1+=1
                elif -15 < j[1] < 0:
                    codigo.insert(cnt1,'bx')
                    cnt1+=1
            elif -30.1631 > j[0] >= -40.5602:
                codigo.insert(cnt1,'AO')
                cnt1+=1                
                if 0 <= j[1] <= 12:
                    codigo.insert(cnt1,'ba')
                    cnt1+=1
                elif 12 < j[1] <= 24:
                    codigo.insert(cnt1,'bb')
                    cnt1+=1
                elif 24 < j[1] <= 36:
                    codigo.insert(cnt1,'bc')
                    cnt1+=1
                elif 36 < j[1] <= 48:
                    codigo.insert(cnt1,'bd')
                    cnt1+=1
                elif 48 < j[1] <= 60:
                    codigo.insert(cnt1,'be')
                    cnt1+=1
                elif 60 < j[1] <= 72:
                    codigo.insert(cnt1,'bf')
                    cnt1+=1
                elif 72 < j[1] <= 84:
                    codigo.insert(cnt1,'bg')
                    cnt1+=1
                elif 84 < j[1] <= 96:
                    codigo.insert(cnt1,'bh')
                    cnt1+=1
                elif 96 < j[1] <= 108:
                    codigo.insert(cnt1,'bi')
                    cnt1+=1
                elif 108 < j[1] <= 120:
                    codigo.insert(cnt1,'bj')
                    cnt1+=1
                elif 120 < j[1] <= 132:
                    codigo.insert(cnt1,'bk')
                    cnt1+=1
                elif 132 < j[1] <= 144:
                    codigo.insert(cnt1,'bl')
                    cnt1+=1
                elif 144 < j[1] <= 156:
                    codigo.insert(cnt1,'bm')
                    cnt1+=1
                elif 156 < j[1] <= 168:
                    codigo.insert(cnt1,'bn')
                    cnt1+=1
                elif 168 < j[1] <= 180:
                    codigo.insert(cnt1,'bo')
                    cnt1+=1
                elif -180 <= j[1] <= -168:
                    codigo.insert(cnt1,'bp')
                    cnt1+=1
                elif -168 < j[1] <= -156:
                    codigo.insert(cnt1,'bq')
                    cnt1+=1
                elif -156 < j[1] <= -144:
                    codigo.insert(cnt1,'br')
                    cnt1+=1
                elif -144 < j[1] <= -132:
                    codigo.insert(cnt1,'bs')
                    cnt1+=1
                elif -132 < j[1] <= -120:
                    codigo.insert(cnt1,'bt')
                    cnt1+=1
                elif -120 < j[1] <= -108:
                    codigo.insert(cnt1,'bu')
                    cnt1+=1
                elif -108 < j[1] <= -96:
                    codigo.insert(cnt1,'bv')
                    cnt1+=1
                elif -96 < j[1] <= -84:
                    codigo.insert(cnt1,'bw')
                    cnt1+=1
                elif -84 < j[1] <= -72:
                    codigo.insert(cnt1,'bx')
                    cnt1+=1
                elif -72 < j[1] <= 60:
                    codigo.insert(cnt1,'by')
                    cnt1+=1
                elif -60 < j[1] <= -48:
                    codigo.insert(cnt1,'bz')
                    cnt1+=1
                elif -48 < j[1] <= -36:
                    codigo.insert(cnt1,'zy')
                    cnt1+=1
                elif -36 < j[1] <= -24:
                    codigo.insert(cnt1,'xw')
                    cnt1+=1
                elif -24 < j[1] <= -12:
                    codigo.insert(cnt1,'vu')
                    cnt1+=1
                elif -12 < j[1] < 0:
                    codigo.insert(cnt1,'ts')
                    cnt1+=1
            elif -20.7738 > j[0] >= -30.1631:
                codigo.insert(cnt1,'AP')
                cnt1+=1                 
                if 0 <= j[1] <= 12:
                    codigo.insert(cnt1,'ba')
                    cnt1+=1
                elif 12 < j[1] <= 24:
                    codigo.insert(cnt1,'bb')
                    cnt1+=1
                elif 24 < j[1] <= 36:
                    codigo.insert(cnt1,'bc')
                    cnt1+=1
                elif 36 < j[1] <= 48:
                    codigo.insert(cnt1,'bd')
                    cnt1+=1
                elif 48 < j[1] <= 60:
                    codigo.insert(cnt1,'be')
                    cnt1+=1
                elif 60 < j[1] <= 72:
                    codigo.insert(cnt1,'bf')
                    cnt1+=1
                elif 72 < j[1] <= 84:
                    codigo.insert(cnt1,'bg')
                    cnt1+=1
                elif 84 < j[1] <= 96:
                    codigo.insert(cnt1,'bh')
                    cnt1+=1
                elif 96 < j[1] <= 108:
                    codigo.insert(cnt1,'bi')
                    cnt1+=1
                elif 108 < j[1] <= 120:
                    codigo.insert(cnt1,'bj')
                    cnt1+=1
                elif 120 < j[1] <= 132:
                    codigo.insert(cnt1,'bk')
                    cnt1+=1
                elif 132 < j[1] <= 144:
                    codigo.insert(cnt1,'bl')
                    cnt1+=1
                elif 144 < j[1] <= 156:
                    codigo.insert(cnt1,'bm')
                    cnt1+=1
                elif 156 < j[1] <= 168:
                    codigo.insert(cnt1,'bn')
                    cnt1+=1
                elif 168 < j[1] <= 180:
                    codigo.insert(cnt1,'bo')
                    cnt1+=1
                elif -180 <= j[1] <= -168:
                    codigo.insert(cnt1,'bp')
                    cnt1+=1
                elif -168 < j[1] <= -156:
                    codigo.insert(cnt1,'bq')
                    cnt1+=1
                elif -156 < j[1] <= -144:
                    codigo.insert(cnt1,'br')
                    cnt1+=1
                elif -144 < j[1] <= -132:
                    codigo.insert(cnt1,'bs')
                    cnt1+=1
                elif -132 < j[1] <= -120:
                    codigo.insert(cnt1,'bt')
                    cnt1+=1
                elif -120 < j[1] <= -108:
                    codigo.insert(cnt1,'bu')
                    cnt1+=1
                elif -108 < j[1] <= -96:
                    codigo.insert(cnt1,'bv')
                    cnt1+=1
                elif -96 < j[1] <= -84:
                    codigo.insert(cnt1,'bw')
                    cnt1+=1
                elif -84 < j[1] <= -72:
                    codigo.insert(cnt1,'bx')
                    cnt1+=1
                elif -72 < j[1] <= 60:
                    codigo.insert(cnt1,'by')
                    cnt1+=1
                elif -60 < j[1] <= -48:
                    codigo.insert(cnt1,'bz')
                    cnt1+=1
                elif -48 < j[1] <= -36:
                    codigo.insert(cnt1,'zy')
                    cnt1+=1
                elif -36 < j[1] <= -24:
                    codigo.insert(cnt1,'xw')
                    cnt1+=1
                elif -24 < j[1] <= -12:
                    codigo.insert(cnt1,'vu')
                    cnt1+=1
                elif -12 < j[1] < 0:
                    codigo.insert(cnt1,'ts')
                    cnt1+=1
            elif -10.2148 > j[0] >= -20.7738:
                codigo.insert(cnt1,'AQ')
                cnt1+=1                
                if 0 <= j[1] <= 10:
                    codigo.insert(cnt1,'ba')
                    cnt1+=1
                elif 10 < j[1] <= 20:
                    codigo.insert(cnt1,'bb')
                    cnt1+=1
                elif 20 < j[1] <= 30:
                    codigo.insert(cnt1,'bc')
                    cnt1+=1
                elif 30 < j[1] <= 40:
                    codigo.insert(cnt1,'bd')
                    cnt1+=1
                elif 40 < j[1] <= 50:
                    codigo.insert(cnt1,'be')
                    cnt1+=1
                elif 50 < j[1] <= 60:
                    codigo.insert(cnt1,'bf')
                    cnt1+=1
                elif 60 < j[1] <= 70:
                    codigo.insert(cnt1,'bg')
                    cnt1+=1
                elif 70 < j[1] <= 80:
                    codigo.insert(cnt1,'bh')
                    cnt1+=1
                elif 80 < j[1] <= 90:
                    codigo.insert(cnt1,'bi')
                    cnt1+=1
                elif 90 < j[1] <= 100:
                    codigo.insert(cnt1,'bj')
                    cnt1+=1
                elif 100 < j[1] <= 110:
                    codigo.insert(cnt1,'bk')
                    cnt1+=1
                elif 110 < j[1] <= 120:
                    codigo.insert(cnt1,'bl')
                    cnt1+=1
                elif 120 < j[1] <= 130:
                    codigo.insert(cnt1,'bm')
                    cnt1+=1
                elif 130 < j[1] <= 140:
                    codigo.insert(cnt1,'bn')
                    cnt1+=1
                elif 140 < j[1] <= 150:
                    codigo.insert(cnt1,'bo')
                    cnt1+=1
                elif 150 < j[1] <= 160:
                    codigo.insert(cnt1,'bp')
                    cnt1+=1
                elif 160 < j[1] <= 170:
                    codigo.insert(cnt1,'bq')
                    cnt1+=1
                elif 170 < j[1] <= 180:
                    codigo.insert(cnt1,'br')
                    cnt1+=1
                elif -180 <= j[1] <= -170:
                    codigo.insert(cnt1,'bs')
                    cnt1+=1
                elif -170 < j[1] <= -160:
                    codigo.insert(cnt1,'bt')
                    cnt1+=1
                elif -150 < j[1] <= -140:
                    codigo.insert(cnt1,'bu')
                    cnt1+=1
                elif -140 < j[1] <= -130:
                    codigo.insert(cnt1,'bv')
                    cnt1+=1
                elif -130 < j[1] <= -120:
                    codigo.insert(cnt1,'bw')
                    cnt1+=1
                elif -120 < j[1] <= -110:
                    codigo.insert(cnt1,'bx')
                    cnt1+=1
                elif -110 < j[1] <= -100:
                    codigo.insert(cnt1,'by')
                    cnt1+=1
                elif -100 < j[1] <= -90:
                    codigo.insert(cnt1,'bz')
                    cnt1+=1
                elif -90 < j[1] <= -80:
                    codigo.insert(cnt1,'zy')
                    cnt1+=1
                elif -80 < j[1] <= -70:
                    codigo.insert(cnt1,'xw')
                    cnt1+=1
                elif -70 < j[1] <= -60:
                    codigo.insert(cnt1,'vu')
                    cnt1+=1
                elif -60 < j[1] <= -50:
                    codigo.insert(cnt1,'ts')
                    cnt1+=1
                elif -50 < j[1] <= -40:
                    codigo.insert(cnt1,'rq')
                    cnt1+=1
                elif -40 < j[1] <= -30:
                    codigo.insert(cnt1,'po')
                    cnt1+=1
                elif -30 < j[1] <= -20:
                    codigo.insert(cnt1,'nm')
                    cnt1+=1
                elif -20 < j[1] <= -10:
                    codigo.insert(cnt1,'lk')
                    cnt1+=1
                elif -10 < j[1] < 0:
                    codigo.insert(cnt1,'ji')
                    cnt1+=1

            elif 0 > j[0] >= -10.2148:
                codigo.insert(cnt1,'AR')
                cnt1+=1
                if 0 <= j[1] <= 10:
                    codigo.insert(cnt1,'ba')
                    cnt1+=1
                elif 10 < j[1] <= 20:
                    codigo.insert(cnt1,'bb')
                    cnt1+=1
                elif 20 < j[1] <= 30:
                    codigo.insert(cnt1,'bc')
                    cnt1+=1
                elif 30 < j[1] <= 40:
                    codigo.insert(cnt1,'bd')
                    cnt1+=1
                elif 40 < j[1] <= 50:
                    codigo.insert(cnt1,'be')
                    cnt1+=1
                elif 50 < j[1] <= 60:
                    codigo.insert(cnt1,'bf')
                    cnt1+=1
                elif 60 < j[1] <= 70:
                    codigo.insert(cnt1,'bg')
                    cnt1+=1
                elif 70 < j[1] <= 80:
                    codigo.insert(cnt1,'bh')
                    cnt1+=1
                elif 80 < j[1] <= 90:
                    codigo.insert(cnt1,'bi')
                    cnt1+=1
                elif 90 < j[1] <= 100:
                    codigo.insert(cnt1,'bj')
                    cnt1+=1
                elif 100 < j[1] <= 110:
                    codigo.insert(cnt1,'bk')
                    cnt1+=1
                elif 110 < j[1] <= 120:
                    codigo.insert(cnt1,'bl')
                    cnt1+=1
                elif 120 < j[1] <= 130:
                    codigo.insert(cnt1,'bm')
                    cnt1+=1
                elif 130 < j[1] <= 140:
                    codigo.insert(cnt1,'bn')
                    cnt1+=1
                elif 140 < j[1] <= 150:
                    codigo.insert(cnt1,'bo')
                    cnt1+=1
                elif 150 < j[1] <= 160:
                    codigo.insert(cnt1,'bp')
                    cnt1+=1
                elif 160 < j[1] <= 170:
                    codigo.insert(cnt1,'bq')
                    cnt1+=1
                elif 170 < j[1] <= 180:
                    codigo.insert(cnt1,'br')
                    cnt1+=1
                elif -180 <= j[1] <= -170:
                    codigo.insert(cnt1,'bs')
                    cnt1+=1
                elif -170 < j[1] <= -160:
                    codigo.insert(cnt1,'bt')
                    cnt1+=1
                elif -150 < j[1] <= -140:
                    codigo.insert(cnt1,'bu')
                    cnt1+=1
                elif -140 < j[1] <= -130:
                    codigo.insert(cnt1,'bv')
                    cnt1+=1
                elif -130 < j[1] <= -120:
                    codigo.insert(cnt1,'bw')
                    cnt1+=1
                elif -120 < j[1] <= -110:
                    codigo.insert(cnt1,'bx')
                    cnt1+=1
                elif -110 < j[1] <= -100:
                    codigo.insert(cnt1,'by')
                    cnt1+=1
                elif -100 < j[1] <= -90:
                    codigo.insert(cnt1,'bz')
                    cnt1+=1
                elif -90 < j[1] <= -80:
                    codigo.insert(cnt1,'zy')
                    cnt1+=1
                elif -80 < j[1] <= -70:
                    codigo.insert(cnt1,'xw')
                    cnt1+=1
                elif -70 < j[1] <= -60:
                    codigo.insert(cnt1,'vu')
                    cnt1+=1
                elif -60 < j[1] <= -50:
                    codigo.insert(cnt1,'ts')
                    cnt1+=1
                elif -50 < j[1] <= -40:
                    codigo.insert(cnt1,'rq')
                    cnt1+=1
                elif -40 < j[1] <= -30:
                    codigo.insert(cnt1,'po')
                    cnt1+=1
                elif -30 < j[1] <= -20:
                    codigo.insert(cnt1,'nm')
                    cnt1+=1
                elif -20 < j[1] <= -10:
                    codigo.insert(cnt1,'lk')
                    cnt1+=1
                elif -10 < j[1] < 0:
                    codigo.insert(cnt1,'ji')
                    cnt1+=1
        codigo.insert(cnt1,'X') #se inserta una X para identificar los frames, por separado
        cnt1+=1                                       
    ss, sss = '', ''    #creando strings 
    for i in codigo:
        ss = i
        sss = sss+ss    #concatena strings de cada frame
    return codigo, sss  #devuelve un string de la secuencia, cada frame concatenado 


def main():

    localA, localB = 'localesA.csv','localesB.csv'   
    posEsfA = conversion(localA) #arreglo de posiciones locales de A en cita y phi
    posEsfB = conversion(localB) #arreglo de posiciones locales de B en cita y phi
    
    cod, s = particion(posEsfA)   #usando funcion particion, para secuencia A 
    file = open('textA.txt', 'w') #guardando string de secuencia A, en textA.txt
    file.write(s)
    file.close()
    
    cod2, s2 = particion(posEsfB) #usando funcion particion, para secuencia B   
    file= open('textB.txt', 'w')  #guardando string de secuencia B, en textB.txt
    file.write(s2)
    file.close()
    
main = main()
