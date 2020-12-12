#Step6.py "4D sequence time alignment" 
#Johnatan Quintero V. johnatan2425@hotmail.com

#Step6: Toma los archivos B_hierarchy.csv, newB_pos.csv y newB_rot.csv
#       Se genera la secuencia C.bvh

import os                #para ejecutar comandos en bash 

def main():

    command = "csv2bvh -o ./C.bvh B_hierarchy.csv newB_pos.csv newB_rot.csv" 
    print("Ejecutando el comando: \n", command)  #command crea el archivo C.bvh
    os.system(command)  #ejecutando comando en terminal

main = main()     


